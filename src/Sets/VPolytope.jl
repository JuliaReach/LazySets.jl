import Base: rand,
             ∈

export VPolytope,
       vertices_list,
       linear_map,
       remove_redundant_vertices,
       tohrep,
       tovrep

"""
    VPolytope{N, VN<:AbstractVector{N}, VT<:AbstractVector{VN}} <: AbstractPolytope{N}

Type that represents a convex polytope in vertex representation.

### Fields

- `vertices` -- list of vertices

### Examples

A polytope in vertex representation can be constructed by passing the list of
vertices. For example, we can build the tetrahedron:

```jldoctest
julia> P = VPolytope([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]);

julia> P.vertices
4-element Vector{Vector{Int64}}:
 [0, 0, 0]
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
```

Alternatively, a `VPolytope` can be constructed passing a matrix of vertices,
where each *column* represents a vertex:

```jldoctest
julia> M = [0 0 0; 1 0 0; 0 1 0; 0 0 1]'
3×4 adjoint(::Matrix{Int64}) with eltype Int64:
 0  1  0  0
 0  0  1  0
 0  0  0  1

julia> P = VPolytope(M);

julia> P.vertices
4-element Vector{Vector{Int64}}:
 [0, 0, 0]
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
```
"""
struct VPolytope{N,VN<:AbstractVector{N},VT<:AbstractVector{VN}} <: AbstractPolytope{N}
    vertices::VT
end

isoperationtype(::Type{<:VPolytope}) = false

# constructor with empty vertices list
VPolytope{N}() where {N} = VPolytope(Vector{Vector{N}}())

# constructor with no vertices of type Float64
VPolytope() = VPolytope{Float64}()

# constructor from rectangular matrix
function VPolytope(vertices_matrix::MT) where {N,MT<:AbstractMatrix{N}}
    vertices = [vertices_matrix[:, j] for j in axes(vertices_matrix, 2)]
    return VPolytope(vertices)
end

"""
    dim(P::VPolytope)

Return the dimension of a polytope in vertex representation.

### Input

- `P` -- polytope in vertex representation

### Output

The ambient dimension of the polytope in vertex representation.
If `P` is empty, the result is ``-1``.

### Examples

```jldoctest
julia> P = VPolytope();

julia> isempty(P.vertices)
true

julia> dim(P)
-1

julia> P = VPolytope([ones(3)]);

julia> P.vertices
1-element Vector{Vector{Float64}}:
 [1.0, 1.0, 1.0]

julia> dim(P) == 3
true
```
"""
function dim(P::VPolytope)
    return isempty(P.vertices) ? -1 : @inbounds length(P.vertices[1])
end

"""
    σ(d::AbstractVector, P::VPolytope)

Return a support vector of a polytope in vertex representation in a given
direction.

### Input

- `d` -- direction
- `P` -- polytope in vertex representation

### Output

A support vector in the given direction.

### Algorithm

A support vector maximizes the support function.
For a polytope, the support function is always maximized in some vertex.
Hence it is sufficient to check all vertices.
"""
function σ(d::AbstractVector, P::VPolytope)
    return _σ_vertices(d, P.vertices)
end

function _σ_vertices(d, vlist)
    # base cases
    m = length(vlist)
    if m == 0
        error("the support vector of an empty polytope is undefined")
    elseif m == 1
        @inbounds return vlist[1]
    end

    # evaluate support function in every vertex
    N = promote_type(eltype(d), eltype(@inbounds vlist[1]))
    max_ρ = N(-Inf)
    max_idx = 0
    for (i, vi) in enumerate(vlist)
        ρ_i = dot(d, vi)
        if ρ_i > max_ρ
            max_ρ = ρ_i
            max_idx = i
        end
    end
    @inbounds return vlist[max_idx]
end

"""
    ρ(d::AbstractVector, P::VPolytope)

Evaluate the support function of a polytope in vertex representation in a given
direction.

### Input

- `d` -- direction
- `P` -- polytope in vertex representation

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, P::VPolytope)
    return _ρ_vertices(d, P.vertices)
end

function _ρ_vertices(d, vlist)
    if isempty(vlist)
        error("the support function of an empty polytope is undefined")
    end
    # evaluate support function in every vertex
    return maximum(v -> dot(d, v), vlist)
end

"""
    ∈(x::AbstractVector{N}, P::VPolytope{N};
      solver=default_lp_solver(N)) where {N}

Check whether a given point is contained in a polytope in vertex representation.

### Input

- `x`      -- point/vector
- `P`      -- polytope in vertex representation
- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

### Output

`true` iff ``x ∈ P``.

### Algorithm

We check, using linear programming, the definition of a convex polytope that a
point is in the set if and only if it is a convex combination of the vertices.

Let ``\\{v_j\\}`` be the ``m`` vertices of `P`.
Then we solve the following ``m``-dimensional linear program.

```math
\\max 0 \\text{ s.t. }
⋀_{i=1}^n ∑_{j=1}^m λ_j v_j[i] = x[i]
∧ ∑_{j=1}^m λ_j = 1
∧ ⋀_{j=1}^m λ_j ≥ 0
```
"""
function ∈(x::AbstractVector{N}, P::VPolytope{N};
           solver=default_lp_solver(N)) where {N}
    vertices = P.vertices
    m = length(vertices)

    # special cases: 0 or 1 vertex
    if m == 0
        return false
    elseif m == 1
        return isapprox(x, @inbounds vertices[1])
    end

    n = length(x)
    @assert n == dim(P) "a vector of length $(length(x)) cannot be " *
                        "contained in a polytope of dimension $(dim(P))"

    A = Matrix{N}(undef, n + 1, m)
    for (j, v_j) in enumerate(vertices)
        # ⋀_i Σ_j λ_j v_j[i] = x[i]
        for i in 1:n
            A[i, j] = v_j[i]
        end
        # Σ_j λ_j = 1
        A[n + 1, j] = one(N)
    end
    b = [x; one(N)]
    lbounds = zero(N)
    ubounds = N(Inf)
    sense = '='
    obj = zeros(N, m)
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    if is_lp_optimal(lp.status)
        return true
    elseif is_lp_infeasible(lp.status)
        return false
    end
    return error("LP returned status $(lp.status) unexpectedly")
end

"""
    rand(::Type{VPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_vertices]::Int=-1)

Create a random polytope in vertex representation.

### Input

- `VPolytope`    -- type for dispatch
- `N`            -- (optional, default: `Float64`) numeric type
- `dim`          -- (optional, default: 2) dimension
- `rng`          -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`         -- (optional, default: `nothing`) seed for reseeding
- `num_vertices` -- (optional, default: `-1`) upper bound on the number of
                    vertices of the polytope (see comment below)

### Output

A random polytope in vertex representation.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of vertices can be controlled with the argument `num_vertices`.
For a negative value we choose a random number in the range `dim:5*dim` (except
if `dim == 1`, in which case we choose in the range `1:2`).
Note that we do not guarantee that the vertices are not redundant.
"""
function rand(::Type{VPolytope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_vertices::Int=-1)
    rng = reseed!(rng, seed)
    if num_vertices < 0
        num_vertices = (dim == 1) ? rand(rng, 1:2) : rand(rng, dim:(5 * dim))
    end
    vertices = [randn(rng, N, dim) for i in 1:num_vertices]
    return VPolytope(vertices)
end

"""
    linear_map(M::AbstractMatrix, P::VPolytope; [apply_convex_hull]::Bool=false)

Concrete linear map of a polytope in vertex representation.

### Input

- `M` -- matrix
- `P` -- polytope in vertex representation
- `apply_convex_hull` -- (optional, default: `false`) flag for applying a convex
                         hull to eliminate redundant vertices

### Output

A polytope in vertex representation.

### Algorithm

The linear map ``M`` is applied to each vertex of the given set ``P``, obtaining
a polytope in vertex representation. The output type is again a `VPolytope`.
"""
function linear_map(M::AbstractMatrix, P::VPolytope;
                    apply_convex_hull::Bool=false)
    @assert dim(P) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(P))"

    return _linear_map_vrep(M, P; apply_convex_hull=apply_convex_hull)
end

@inline function _linear_map_vrep(M::AbstractMatrix, P::VPolytope,
                                  algo::LinearMapVRep=LinearMapVRep(nothing);  # ignored
                                  apply_convex_hull::Bool=false)
    vlist = broadcast(v -> M * v, P.vertices)
    if apply_convex_hull
        convex_hull!(vlist)
    end
    return VPolytope(vlist)
end

"""
    translate(P::VPolytope, v::AbstractVector)

Translate (i.e., shift) a polytope in vertex representation by a given vector.

### Input

- `P` -- polytope in vertex representation
- `v` -- translation vector

### Output

A translated polytope in vertex representation.

### Notes

See [`translate!(::VPolytope, ::AbstractVector)`](@ref) for the in-place version.
"""
function translate(P::VPolytope, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

"""
    translate!(P::VPolytope, v::AbstractVector)

Translate (i.e., shift) a polytope in vertex representation by a given vector,
in-place.

### Input

- `P` -- polytope in vertex representation
- `v` -- translation vector

### Output

The polytope `P` translated by `v`.

### Notes

See [`translate(::VPolytope, ::AbstractVector)`](@ref) for the out-of-place
version.

### Algorithm

We add the vector to each vertex of the polytope.
"""
function translate!(P::VPolytope, v::AbstractVector)
    if isempty(P.vertices)
        return P
    end

    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    for x in P.vertices
        x .+= v
    end
    return P
end

"""
    vertices_list(P::VPolytope)

Return the list of vertices of a polytope in vertex representation.

### Input

- `P` -- polytope in vertex representation

### Output

The list of vertices.
"""
function vertices_list(P::VPolytope)
    return P.vertices
end

"""
    constraints_list(P::VPolytope)

Return a list of constraints defining a polytope in vertex representation.

### Input

- `P` -- polytope in vertex representation

### Output

A list of constraints of the polytope.

### Algorithm

We use `tohrep` to compute the constraint representation of `P`.
"""
function constraints_list(P::VPolytope)
    n = dim(P)
    if n == 1
        return constraints_list(convert(Interval, P))
    elseif n == 2
        return constraints_list(convert(VPolygon, P))
    else
        return constraints_list(tohrep(P))
    end
end

"""
    remove_redundant_vertices(P::VPolytope{N};
                              [backend]=nothing,
                              [solver]=nothing) where {N}

Return the polytope obtained by removing the redundant vertices of the given
polytope in vertex representation.

### Input

- `P`       -- polytope in vertex representation
- `backend` -- (optional, default: `nothing`) the backend for polyhedral
               computations; see `default_polyhedra_backend(P)` or
               [Polyhedra's documentation](https://juliapolyhedra.github.io/)
               for further information
- `solver`  -- (optional, default: `nothing`) the linear programming
               solver used in the backend, if needed; see
               `default_lp_solver_polyhedra(N)`

### Output

A new polytope such that its vertices are the convex hull of the given polytope.

### Notes

The optimization problem associated to removing redundant vertices is handled
by `Polyhedra`.
If the polyhedral computations backend requires an LP solver, but it has not
been specified, we use `default_lp_solver_polyhedra(N)` to define such solver.
Otherwise, the redundancy-removal function of the polyhedral backend is used.
"""
function remove_redundant_vertices(P::VPolytope{N};
                                   backend=nothing,
                                   solver=nothing) where {N}
    require(@__MODULE__, :Polyhedra; fun_name="remove_redundant_vertices")
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end
    Q = polyhedron(P; backend=backend)
    if Polyhedra.supportssolver(typeof(Q))
        if isnothing(solver)
            # presolver prints warnings about infeasible solutions (#3226)
            solver = default_lp_solver_polyhedra(N; presolve=false)
        end
        vQ = Polyhedra.vrep(Q)
        Polyhedra.setvrep!(Q, Polyhedra.removevredundancy(vQ, solver))
    else
        removevredundancy!(Q; ztol=_ztol(N))
    end
    return VPolytope(Q)
end

"""
    tohrep(P::VPolytope{N};
           [backend]=default_polyhedra_backend(P)) where {N}

Transform a polytope in vertex representation to a polytope in constraint
representation.

### Input

- `P`       -- polytope in vertex representation
- `backend` -- (optional, default: `default_polyhedra_backend(P)`) the
               backend for polyhedral computations; see [Polyhedra's
               documentation](https://juliapolyhedra.github.io/) for further
               information

### Output

A `HPolytope` as the constraint representation of `P`.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.
"""
function tohrep(P::VPolytope{N};
                backend=default_polyhedra_backend(P)) where {N}
    @assert !isempty(P.vertices) "cannot convert an empty polytope in vertex " *
                                 "representation to constraint representation"
    require(@__MODULE__, :Polyhedra; fun_name="tohrep")
    return convert(HPolytope, polyhedron(P; backend=backend))
end

"""
    tovrep(P::VPolytope)

Return a vertex representation of the given polytope in vertex representation
(no-op).

### Input

- `P` -- polytope in vertex representation

### Output

The same polytope instance.
"""
function tovrep(P::VPolytope)
    return P
end

function load_polyhedra_vpolytope() # function to be loaded by Requires
    return quote
        # see the interface file init_Polyhedra.jl for the imports

        # VPolytope from a VRep
        function VPolytope(P::VRep{N}) where {N}
            vertices = collect(Polyhedra.points(P))
            return VPolytope(vertices)
        end

        """
            polyhedron(P::VPolytope;
                       [backend]=default_polyhedra_backend(P),
                       [relative_dimension]=nothing)

        Return a `VRep` polyhedron from `Polyhedra.jl` given a polytope in vertex
        representation.

        ### Input

        - `P`       -- polytope in vertex representation
        - `backend` -- (optional, default: `default_polyhedra_backend(P)`) the
                       backend for polyhedral computations; see [Polyhedra's
                       documentation](https://juliapolyhedra.github.io/) for further
                       information
        - `relative_dimension` -- (default, optional: `nothing`) an integer representing
                                  the (relative) dimension of the polytope; this
                                  argument is mandatory if the polytope is empty

        ### Output

        A `VRep` polyhedron.

        ### Notes

        The *relative dimension* (or just *dimension*) refers to the dimension of the
        set relative to itself, independently of the ambient dimension. For example, a
        point has (relative) dimension zero, and a line segment has (relative) dimension
        one.

        In this library, `LazySets.dim` always returns the ambient dimension of the set,
        such that a line segment in two dimensions has dimension two. However,
        `Polyhedra.dim` will assign a dimension equal to one to a line segment
        because it uses a different convention.
        """
        function polyhedron(P::VPolytope;
                            backend=default_polyhedra_backend(P),
                            relative_dimension=nothing)
            if isempty(P)
                if isnothing(relative_dimension)
                    error("the conversion to a `Polyhedra.polyhedron` requires the " *
                          "(relative) dimension of the `VPolytope` to be known, but it " *
                          "cannot be inferred from an empty set; use the keyword " *
                          "argument `relative_dimension`")
                end
                return polyhedron(Polyhedra.vrep(P.vertices; d=relative_dimension),
                                  backend)
            end
            return polyhedron(Polyhedra.vrep(P.vertices), backend)
        end
    end
end  # quote / load_polyhedra_vpolytope()

function project(P::VPolytope, block::AbstractVector{Int}; kwargs...)
    if isempty(P.vertices)
        return P
    end

    m = length(block)
    if m == 1
        l, h = extrema(P, block[1])
        return Interval(l, h)
    end

    M = projection_matrix(block, dim(P), eltype(P.vertices))
    πvertices = broadcast(v -> M * v, P.vertices)

    if m == 2
        return VPolygon(πvertices; apply_convex_hull=true)
    else
        return VPolytope(πvertices)
    end
end

function permute(P::VPolytope, p::AbstractVector{Int})
    vlist = [v[p] for v in P.vertices]
    return VPolytope(vlist)
end

"""
    reflect(P::VPolytope)

Concrete reflection of a polytope in vertex representation `P`, resulting in the
reflected set `-P`.

### Input

- `P` -- polytope in vertex representation

### Output

The `VPolytope` representing `-P`.
"""
function reflect(P::VPolytope)
    return VPolytope(-P.vertices)
end

function scale!(α::Real, P::VPolytope)
    P.vertices .*= α
    return P
end
