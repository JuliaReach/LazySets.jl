using MathProgBase, GLPKMathProgInterface

import Base: rand,
             ∈

export VPolytope,
       vertices_list,
       convex_hull,
       cartesian_product,
       linear_map,
       remove_redundant_vertices,
       minkowski_sum

"""
    VPolytope{N<:Real} <: AbstractPolytope{N}

Type that represents a convex polytope in V-representation.

### Fields

- `vertices` -- the list of vertices
"""
struct VPolytope{N<:Real} <: AbstractPolytope{N}
    vertices::Vector{Vector{N}}
end

# constructor for a VPolytope with empty vertices list
VPolytope{N}() where {N<:Real} = VPolytope{N}(Vector{Vector{N}}())

# constructor for a VPolytope with no vertices of type Float64
VPolytope() = VPolytope{Float64}()

# constructor from a polygon in V-representation
function VPolytope(P::VPolygon, share::Bool=false)
    v = vertices_list(P)
    if !share
        v = copy(v)
    end
    return VPolytope(v)
end

# --- LazySet interface functions ---

"""
    dim(P::VPolytope)::Int

Return the dimension of a polytope in V-representation.

### Input

- `P`  -- polytope in V-representation

### Output

The ambient dimension of the polytope in V-representation.
If it is empty, the result is ``-1``.

### Examples

```jldoctest
julia> v = VPolytope();

julia> dim(v) > 0
false

julia> v = VPolytope([ones(3)])
VPolytope{Float64}(Array{Float64,1}[[1.0, 1.0, 1.0]])

julia> dim(v) == 3
true

```
"""
function dim(P::VPolytope)::Int
    return length(P.vertices) == 0 ? -1 : length(P.vertices[1])
end

"""
    σ(d::AbstractVector{N}, P::VPolytope{N}) where {N<:Real}

Return the support vector of a polytope in V-representation in a given
direction.

### Input

- `d` -- direction
- `P` -- polytope in V-representation

### Output

The support vector in the given direction.

### Algorithm

A support vector maximizes the support function.
For a polytope, the support function is always maximized in some vertex.
Hence it is sufficient to check all vertices.
"""
function σ(d::AbstractVector{N}, P::VPolytope{N}) where {N<:Real}
    # base cases
    m = length(P.vertices)
    if m == 0
        error("the support function for an empty polytope is undefined")
    elseif m == 1
        return P.vertices[1]
    end

    # evaluate support function in every vertex
    max_ρ = N(-Inf)
    max_idx = 0
    for (i, vi) in enumerate(P.vertices)
        ρ_i = dot(d, vi)
        if ρ_i > max_ρ
            max_ρ = ρ_i
            max_idx = i
        end
    end
    return P.vertices[max_idx]
end

"""
    ∈(x::AbstractVector{N}, P::VPolytope{N};
      solver=default_lp_solver(N))::Bool where {N<:Real}

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
\\max 0 \\text{s.t.}
\bigwedge_{i=1}^n \\sum_{j=1}^m λ_j v_j[i] = x[i]
\\sum_{j=1}^m λ_j = 1
\bigwedge_{j=1}^m λ_j ≥ 0
```
"""
function ∈(x::AbstractVector{N}, P::VPolytope{N};
           solver=default_lp_solver(N))::Bool where {N<:Real}
    vertices = P.vertices
    m = length(vertices)

    # special cases: 0 or 1 vertex
    if m == 0
        return false
    elseif m == 1
        return x == vertices[1]
    end

    n = length(x)
    @assert n == dim(P) "a vector of length $(length(x)) cannot be " *
        "contained in a polytope of dimension $(dim(P))"

    A = Matrix{N}(undef, n+1, m)
    for j in 1:m
        v_j = vertices[j]
        # ⋀_i Σ_j λ_j v_j[i] = x[i]
        for i in 1:n
            A[i, j] = v_j[i]
        end
        # Σ_j λ_j = 1
        A[n+1, j] = one(N)
    end
    b = [x; one(N)]

    lbounds = zeros(N, m)
    ubounds = Inf
    sense = [fill('=', n); '<']
    obj = zeros(N, m)

    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    if lp.status == :Optimal
        return true
    elseif lp.status == :Infeasible
        return false
    end
    @assert false "LP returned status $(lp.status) unexpectedly"
end

"""
    rand(::Type{VPolytope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_vertices]::Int=-1)::VPolytope{N}

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
              seed::Union{Int, Nothing}=nothing,
              num_vertices::Int=-1
             )::VPolytope{N}
    rng = reseed(rng, seed)
    if num_vertices < 0
        num_vertices = (dim == 1) ? rand(1:2) : rand(dim:5*dim)
    end
    vertices = [randn(rng, N, dim) for i in 1:num_vertices]
    return VPolytope(vertices)
end

"""
    linear_map(M::AbstractMatrix{N}, P::VPolytope{N}) where {N<:Real}

Concrete linear map of a polytope in vertex representation.

### Input

- `M` -- matrix
- `P` -- polytope in vertex representation

### Output

A polytope in vertex representation.

### Algorithm

The linear map ``M`` is applied to each vertex of the given set ``P``, obtaining
a polytope in V-representation. The output type is again a `VPolytope`.
"""
function linear_map(M::AbstractMatrix{N}, P::VPolytope{N}) where {N<:Real}
    @assert dim(P) == size(M, 2) "a linear map of size $(size(M)) cannot be applied to a set of dimension $(dim(P))"
    return _linear_map_vrep(M, P)
end

@inline function _linear_map_vrep(M::AbstractMatrix{N}, P::VPolytope{N}) where {N<:Real}
    return broadcast(v -> M * v, vertices_list(P)) |> VPolytope{N}
end

"""
    translate(P::VPolytope{N}, v::AbstractVector{N}) where {N<:Real}

Translate (i.e., shift) a polytope in vertex representation by a given vector.

### Input

- `P` -- polytope in vertex representation
- `v` -- translation vector

### Output

A translated polytope in vertex representation.

### Algorithm

We add the vector to each vertex of the polytope.
"""
function translate(P::VPolytope{N}, v::AbstractVector{N}) where {N<:Real}
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return VPolytope([x + v for x in vertices_list(P)])
end

# --- AbstractPolytope interface functions ---

"""
    vertices_list(P::VPolytope{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polytope in V-representation.

### Input

- `P` -- polytope in vertex representation

### Output

List of vertices.
"""
function vertices_list(P::VPolytope{N})::Vector{Vector{N}} where {N<:Real}
    return P.vertices
end

"""
    constraints_list(P::VPolytope{N}) where {N<:Real}

Return the list of constraints defining a polytope in V-representation.

### Input

- `P` -- polytope in V-representation

### Output

The list of constraints of the polytope.

### Algorithm

First the H-representation of ``P`` is computed, then its list of constraints
is returned. 
"""
function constraints_list(P::VPolytope{N}) where {N<:Real}
    return constraints_list(tohrep(P))
end

"""
    convex_hull(P1::VPolytope{N}, P2::VPolytope{N};
                [backend]=default_polyhedra_backend(P1, N)) where {N<:Real}

Compute the convex hull of the set union of two polytopes in V-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `nothing`) the polyhedral
                  computations backend

### Output

The `VPolytope` obtained by the concrete convex hull of `P1` and `P2`.

### Notes

This function takes the union of the vertices of each polytope and then relies
on a concrete convex hull algorithm. For low dimensions, a specialized implementation
for polygons is used. For higher dimensions, `convex_hull` relies on the polyhedral
computations backend that can be specified using the `backend` keyword argument.

For performance reasons, it is suggested to use the `CDDLib.Library()` backend
for the `convex_hull`.
"""
function convex_hull(P1::VPolytope{N}, P2::VPolytope{N};
                     backend=nothing) where {N<:Real}
    vunion = [P1.vertices; P2.vertices]
    convex_hull!(vunion; backend=backend)
    return VPolytope(vunion)
end

"""
    cartesian_product(P1::VPolytope{N}, P2::VPolytope{N};
                      [backend]=default_polyhedra_backend(P1, N)) where {N}

Compute the Cartesian product of two polytopes in V-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `default_polyhedra_backend(P1, N)`) the polyhedral
                  computations backend, see
                  [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `VPolytope` obtained by the concrete Cartesian product of `P1` and `P2`.
"""
function cartesian_product(P1::VPolytope{N}, P2::VPolytope{N};
                           backend=default_polyhedra_backend(P1, N)) where {N}
    require(:Polyhedra; fun_name="cartesian_product")
    Pcp = Polyhedra.vcartesianproduct(polyhedron(P1; backend=backend),
                                      polyhedron(P2; backend=backend))
    return VPolytope(Pcp)
end

"""
    remove_redundant_vertices(P::VPolytope{N};
                              [backend]=nothing,
                              [solver]=nothing)::VPolytope{N} where {N<:Real}

Return the polytope obtained by removing the redundant vertices of the given polytope.

### Input

- `P`       -- polytope in vertex representation
- `backend` -- (optional, default: `nothing`) the backend for polyhedral
               computations; see `default_polyhedra_backend(P, N)` or
               [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
               for further information on the available backends
- `solver`  -- (optional, default: `nothing`) the linear programming
               solver used in the backend, if needed; see
               `default_lp_solver_polyhedra(N)`

### Output

A new polytope such that its vertices are the convex hull of the given polytope.

### Notes

The optimization problem associated to removing redundant vertices is handled
by `Polyhedra`.
If the polyhedral computations backend requires an LP solver but it has not been
set, we use `default_lp_solver_polyhedra(N)` to define such solver.
Otherwise, the redundancy removal function of the polyhedral backend is used.
"""
function remove_redundant_vertices(P::VPolytope{N};
                                   backend=nothing,
                                   solver=nothing)::VPolytope{N} where {N<:Real}
    require(:Polyhedra; fun_name="remove_redundant_vertices")
    if backend == nothing
        backend = default_polyhedra_backend(P, N)
    end
    Q = polyhedron(P; backend=backend)
    if Polyhedra.supportssolver(typeof(Q))
        if solver == nothing
            solver = default_lp_solver_polyhedra(N)
        end
        vQ = Polyhedra.vrep(Q)
        Polyhedra.setvrep!(Q, Polyhedra.removevredundancy(vQ, solver))
    else
        Polyhedra.removevredundancy!(Q)
    end
    return VPolytope(Q)
end

"""
    tohrep(P::VPolytope{N};
           [backend]=default_polyhedra_backend(P, N)) where {N<:Real}

Transform a polytope in V-representation to a polytope in H-representation.

### Input

- `P`       -- polytope in vertex representation
- `backend` -- (optional, default: `default_polyhedra_backend(P, N)`) the
               backend for polyhedral computations

### Output

The `HPolytope` which is the constraint representation of the given polytope
in vertex representation.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.
For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1).
"""
function tohrep(P::VPolytope{N};
                backend=default_polyhedra_backend(P, N)) where {N<:Real}
    vl = P.vertices
    if isempty(vl)
        return EmptySet{N}()
    end
    require(:Polyhedra; fun_name="tohrep")
    return HPolytope(polyhedron(P; backend=backend))
end

"""
    tovrep(P::VPolytope)

Return a vertex representation of the given polytope in vertex
representation (no-op).

### Input

- `P` -- polytope in vertex representation

### Output

The same polytope instance.
"""
function tovrep(P::VPolytope)
    return P
end

"""
    minkowski_sum(P1::VPolytope{N}, P2::VPolytope{N};
                  [apply_convex_hull]=true,
                  [backend]=nothing,
                  [solver]=nothing) where {N<:Real}

Compute the Minkowski sum between two polytopes in vertex representation.

### Input

- `P1`                -- polytope
- `P2`                -- another polytope
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         pairwise sums using a convex hull algorithm 
- `backend`           -- (optional, default: `nothing`) the backend for
                         polyhedral computations used to post-process with a
                         convex hull; see `default_polyhedra_backend(P1, N)`
- `solver`            -- (optional, default: `nothing`) the backend used to
                         solve the linear program; see
                         `default_lp_solver_polyhedra(N)`

### Output

A new polytope in vertex representation whose vertices are the convex hull of
the sum of all possible sums of vertices of `P1` and `P2`.
"""
function minkowski_sum(P1::VPolytope{N}, P2::VPolytope{N};
                       apply_convex_hull::Bool=true,
                       backend=nothing,
                       solver=nothing) where {N<:Real}

    @assert dim(P1) == dim(P2) "cannot compute the Minkowski sum between a polyotope " *
        "of dimension $(dim(P1)) and a polytope of dimension $((dim(P2)))"

    vlist1 = _vertices_list(P1, backend)
    vlist2 = _vertices_list(P2, backend)
    n, m = length(vlist1), length(vlist2)
    Vout = Vector{Vector{N}}()
    sizehint!(Vout, n * m)
    for vi in vlist1
        for vj in vlist2
            push!(Vout, vi + vj)
        end
    end
    if apply_convex_hull
        if backend == nothing
            backend = default_polyhedra_backend(P1, N)
            solver = default_lp_solver_polyhedra(N)
        end
        convex_hull!(Vout, backend=backend, solver=solver)
    end
    return VPolytope(Vout)
end

# ==========================================
# Lower level methods that use Polyhedra.jl
# ==========================================

function load_polyhedra_vpolytope() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

export vertices_list,
       tohrep,
       tovrep

# VPolytope from a VRep
function VPolytope(P::VRep{N}) where {N}
    vertices = Vector{Vector{N}}()
    for vi in Polyhedra.points(P)
        push!(vertices, vi)
    end
    return VPolytope(vertices)
end

"""
    polyhedron(P::VPolytope{N};
               [backend]=default_polyhedra_backend(P, N)) where {N<:Real}

Return an `VRep` polyhedron from `Polyhedra.jl` given a polytope in V-representation.

### Input

- `P`       -- polytope
- `backend` -- (optional, default: `default_polyhedra_backend(P, N)`) the polyhedral
               computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
               for further information
- `relative_dimension` -- (default, optional: `nothing`) an integer representing
                          the (relative) dimension of the polytope; this argument
                          is mandatory if the polytope is empty

### Output

A `VRep` polyhedron.

### Notes

The *relative dimension* (or just *dimension*) refers to the dimension of the set
relative to itself, independently of the ambient dimension. For example, a point
has (relative) dimension zero, and a line segment has (relative) dimension one.

In this library, `LazySets.dim` always returns the ambient dimension of the set,
such that a line segment in two dimensions has dimension two. However,
`Polyhedra.dim` will assign a dimension equal to one to a line segment
because it uses a different convention.
"""
function polyhedron(P::VPolytope{N};
                    backend=default_polyhedra_backend(P, N),
                    relative_dimension=nothing) where {N<:Real}
    if isempty(P)
        if relative_dimension == nothing
            error("the conversion to a `Polyhedra.polyhedron` requires the (relative) dimension " *
                  "of the `VPolytope` to be known, but it cannot be inferred from an empty set; " *
                  "try passing it in the keyword argument `relative_dimension`")
        end
        return polyhedron(Polyhedra.vrep(P.vertices, d=relative_dimension), backend)
    end
    return polyhedron(Polyhedra.vrep(P.vertices), backend)
end

end # quote
end # function load_polyhedra_vpolytope()
