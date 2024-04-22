import Base: rand,
             convert

export HPolyhedron,
       dim, σ, ∈,
       addconstraint!,
       constraints_list,
       tohrep, tovrep,
       convex_hull,
       cartesian_product,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       constrained_dimensions,
       is_hyperplanar

"""
    HPolyhedron{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a convex polyhedron in constraint representation, that is,
a finite intersection of half-spaces,
```math
P = ⋂_{i = 1}^m H_i,
```
where each ``H_i = \\{x ∈ ℝ^n : a_i^T x ≤ b_i \\}`` is a
half-space, ``a_i ∈ ℝ^n`` is the normal vector of the ``i``-th
half-space and ``b_i`` is the displacement. The set ``P`` may or may not be
bounded (see also [`HPolytope`](@ref), which assumes boundedness).

### Fields

- `constraints` -- vector of linear constraints
"""
struct HPolyhedron{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    constraints::Vector{HalfSpace{N,VN}}

    function HPolyhedron(constraints::Vector{HalfSpace{N,VN}}) where {N,VN<:AbstractVector{N}}
        return new{N,VN}(constraints)
    end
end

isoperationtype(::Type{<:HPolyhedron}) = false

# constructor with no constraints
function HPolyhedron{N,VN}() where {N,VN<:AbstractVector{N}}
    return HPolyhedron(Vector{HalfSpace{N,VN}}())
end

# constructor with no constraints, given only the numeric type
function HPolyhedron{N}() where {N}
    return HPolyhedron(Vector{HalfSpace{N,Vector{N}}}())
end

# constructor without explicit numeric type, defaults to Float64
function HPolyhedron()
    return HPolyhedron{Float64}()
end

# constructor with constraints of mixed type
function HPolyhedron(constraints::Vector{<:HalfSpace})
    return HPolyhedron(_normal_Vector(constraints))
end

# constructor from a simple constraint representation
function HPolyhedron(A::AbstractMatrix, b::AbstractVector)
    return HPolyhedron(constraints_list(A, b))
end

# convenience union type
const HPoly{N} = Union{HPolytope{N},HPolyhedron{N}}

"""
    dim(P::HPoly)

Return the dimension of a polyhedron in constraint representation.

### Input

- `P`  -- polyhedron in constraint representation

### Output

The ambient dimension of the polyhedron in constraint representation.
If it has no constraints, the result is ``-1``.
"""
function dim(P::HPoly)
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end

"""
    ρ(d::AbstractVector{M}, P::HPoly{N};
      solver=default_lp_solver(M, N)) where {M, N}

Evaluate the support function of a polyhedron in constraint representation in a
given direction.

### Input

- `d`      -- direction
- `P`      -- polyhedron in constraint representation
- `solver` -- (optional, default: `default_lp_solver(M, N)`) the backend used to
              solve the linear program

### Output

The evaluation of the support function for the polyhedron.
If a polytope is unbounded in the given direction, we throw an error.
If a polyhedron is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{M}, P::HPoly{N};
           solver=default_lp_solver(M, N)) where {M,N}
    lp, unbounded = σ_helper(d, P, solver)
    if unbounded
        if P isa HPolytope
            error("the support function in direction $(d) is undefined " *
                  "because the polytope is unbounded")
        end
        return N(Inf)
    end
    return dot(d, lp.sol)
end

"""
    σ(d::AbstractVector{M}, P::HPoly{N};
      solver=default_lp_solver(M, N) where {M, N}

Return a support vector of a polyhedron in constraint representation in a given
direction.

### Input

- `d`      -- direction
- `P`      -- polyhedron in constraint representation
- `solver` -- (optional, default: `default_lp_solver(M, N)`) the backend used to
              solve the linear program

### Output

The support vector in the given direction.
If a polytope is unbounded in the given direction, we throw an error.
If a polyhedron is unbounded in the given direction, the result contains `±Inf`
entries.
"""
function σ(d::AbstractVector{M}, P::HPoly{N};
           solver=default_lp_solver(M, N)) where {M,N}
    lp, unbounded = σ_helper(d, P, solver)
    if unbounded
        if P isa HPolytope
            error("the support vector in direction $(d) is undefined because " *
                  "the polytope is unbounded")
        end
        return _σ_unbounded_lp(d, P, lp)
    else
        return lp.sol
    end
end

# construct the solution from the solver's ray result
function _σ_unbounded_lp(d, P::HPoly{N}, lp) where {N}
    if isnothing(lp)
        ray = d
    elseif has_lp_infeasibility_ray(lp.model)
        ray = lp.sol  # infeasibility ray is stored as the solution
    else
        error("LP solver did not return an infeasibility ray")
    end

    res = Vector{N}(undef, length(ray))
    e = isempty(P.constraints) ? zeros(N, length(ray)) : an_element(P)
    @inbounds for i in eachindex(ray)
        if isapproxzero(ray[i])
            res[i] = e[i]
        elseif ray[i] > zero(N)
            res[i] = N(Inf)
        else
            res[i] = N(-Inf)
        end
    end
    return res
end

function σ_helper(d::AbstractVector, P::HPoly, solver)
    # represent c = -d as a Vector since GLPK does not accept sparse vectors
    # (see #1011)
    c = to_negative_vector(d)

    (A, b) = tosimplehrep(P)
    if length(b) == 0
        unbounded = true
        lp = nothing
    else
        sense = '<'
        l = -Inf
        u = Inf
        lp = linprog(c, A, sense, b, l, u, solver)
        if is_lp_infeasible(lp.status; strict=true)
            throw(ArgumentError("the support vector is undefined because " *
                                "the polyhedron is empty"))
        elseif is_lp_unbounded(lp.status)
            unbounded = true
        elseif is_lp_optimal(lp.status)
            unbounded = false
        else
            error("got unknown LP status $(lp.status)")
        end
    end
    return (lp, unbounded)
end

"""
    rand(::Type{HPolyhedron}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random polyhedron.

### Input

- `HPolyhedron` -- type for dispatch
- `N`           -- (optional, default: `Float64`) numeric type
- `dim`         -- (optional, default: 2) dimension (is ignored)
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A random polyhedron.

### Algorithm

We first create a random polytope and then for each constraint randomly (50%)
decide whether to include it.
"""
function rand(::Type{HPolyhedron};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    P = rand(HPolytope; N=N, dim=dim, rng=rng)
    constraints_P = constraints_list(P)
    constraints_Q = Vector{eltype(constraints_P)}()
    for i in eachindex(constraints_P)
        if rand(rng, Bool)
            push!(constraints_Q, constraints_P[i])
        end
    end
    return HPolyhedron(constraints_Q)
end

"""
    addconstraint!(P::HPoly, constraint::HalfSpace)

Add a linear constraint to a polyhedron in constraint representation.

### Input

- `P`          -- polyhedron in constraint representation
- `constraint` -- linear constraint to add

### Notes

It is left to the user to guarantee that the dimension of all linear constraints
is the same.
"""
function addconstraint!(P::HPoly, constraint::HalfSpace)
    push!(P.constraints, constraint)
    return nothing
end

"""
    constraints_list(P::HPoly)

Return the list of constraints defining a polyhedron in constraint
representation.

### Input

- `P` -- polyhedron in constraint representation

### Output

The list of constraints of the polyhedron.
"""
function constraints_list(P::HPoly)
    return P.constraints
end

"""
    tohrep(P::HPoly)

Return a constraint representation of the given polyhedron in constraint
representation (no-op).

### Input

- `P` -- polyhedron in constraint representation

### Output

The same polyhedron instance.
"""
function tohrep(P::HPoly)
    return P
end

"""
    normalize(P::HPoly{N}, p::Real=N(2)) where {N}

Normalize a polyhedron in constraint representation.

### Input

- `P` -- polyhedron in constraint representation
- `p` -- (optional, default: `2`) norm

### Output

A new polyhedron in constraint representation whose normal directions ``a_i``
are normalized, i.e., such that ``‖a_i‖_p = 1`` holds.
"""
function normalize(P::HPoly{N}, p::Real=N(2)) where {N}
    constraints = [normalize(hs, p) for hs in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end

"""
    remove_redundant_constraints(P::HPoly{N}; [backend]=nothing) where {N}

Remove the redundant constraints in a polyhedron in constraint representation.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `nothing`) the backend used to solve the
               linear program

### Output

A polyhedron equivalent to `P` but with no redundant constraints, or an empty
set if `P` is detected to be empty, which may happen if the constraints are
infeasible.

### Notes

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:HalfSpace})` for details.
"""
function remove_redundant_constraints(P::HPoly; backend=nothing)
    Pred = copy(P)
    if remove_redundant_constraints!(Pred; backend=backend)
        return Pred
    else # the polyhedron P is empty
        N = eltype(P)
        return EmptySet{N}(dim(P))
    end
end

"""
    remove_redundant_constraints!(P::HPoly{N}; [backend]=nothing) where {N}

Remove the redundant constraints of a polyhedron in constraint representation;
the polyhedron is updated in-place.

### Input

- `P`       -- polyhedron
- `backend` -- (optional, default: `nothing`) the backend used to solve the
               linear program

### Output

`true` if the method was successful and the polyhedron `P` is modified by
removing its redundant constraints, and `false` if `P` is detected to be empty,
which may happen if the constraints are infeasible.

### Notes

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:HalfSpace})` for details.
"""
function remove_redundant_constraints!(P::HPoly; backend=nothing)
    return remove_redundant_constraints!(P.constraints; backend=backend)
end

"""
    translate(P::HPoly, v::AbstractVector; [share]::Bool=false)

Translate (i.e., shift) a polyhedron in constraint representation by a given
vector.

### Input

- `P`     -- polyhedron in constraint representation
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated polyhedron in constraint representation.

### Notes

The normal vectors of the constraints (vector `a` in `a⋅x ≤ b`) are shared with
the original constraints if `share == true`.

### Algorithm

We translate every constraint.
"""
function translate(P::HPoly, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    T = basetype(P)
    return T(constraints)
end

"""
    tovrep(P::HPoly; [backend]=default_polyhedra_backend(P))

Transform a polytope in constraint representation to a polytope in vertex
representation.

### Input

- `P`       -- polytope in constraint representation
- `backend` -- (optional, default: `default_polyhedra_backend(P)`) the
               backend for polyhedral computations

### Output

A `VPolytope` which is a vertex representation of the given polytope in
constraint representation.

### Notes

The conversion may not preserve the numeric type (e.g., with `N == Float32`)
depending on the backend.
For further information on the supported backends see [Polyhedra's
documentation](https://juliapolyhedra.github.io/).
"""
function tovrep(P::HPoly; backend=default_polyhedra_backend(P))
    require(@__MODULE__, :Polyhedra; fun_name="tovrep")
    P = polyhedron(P; backend=backend)
    return VPolytope(P)
end

# this method is required mainly for HPolytope (because the fallback for
# AbstractPolytope is incorrect with no constraints)
#
# the method also treats a corner case for problems with Rationals in LP solver
function isempty(P::HPoly,
                 witness::Bool=false;
                 use_polyhedra_interface::Bool=false,
                 solver=nothing,
                 backend=nothing)
    if length(constraints_list(P)) < 2
        return witness ? (false, an_element(P)) : false
    end
    return _isempty_polyhedron(P, witness;
                               use_polyhedra_interface=use_polyhedra_interface,
                               solver=solver, backend=backend)
end

function load_polyhedra_hpolyhedron() # function to be loaded by Requires
    return quote
        # see the interface file init_Polyhedra.jl for the imports

        """
             convert(::Type{HPolyhedron}, P::HRep{N}) where {N}

        Convert an `HRep` polyhedron from `Polyhedra.jl` to a polyhedron in constraint
        representation .

        ### Input

        - `HPolyhedron` -- target type
        - `P`           -- `HRep` polyhedron

        ### Output

        An `HPolyhedron`.
        """
        function convert(::Type{HPolyhedron}, P::HRep{N}) where {N}
            VN = Polyhedra.hvectortype(P)
            constraints = Vector{HalfSpace{N,VN}}()
            for hi in Polyhedra.allhalfspaces(P)
                a, b = hi.a, hi.β
                if isapproxzero(norm(a))
                    continue
                end
                push!(constraints, HalfSpace(a, b))
            end
            return HPolyhedron(constraints)
        end

        # convenience conversion method
        function HPolyhedron(P::HRep{N}) where {N}
            return convert(HPolyhedron, P)
        end
    end
end  # quote / load_polyhedra_hpolyhedron()

function is_hyperplanar(P::HPolyhedron)
    clist = P.constraints
    m = length(clist)

    # check that the number of constraints is fine
    if m > 2
        # try to remove redundant constraints
        clist = remove_redundant_constraints(clist)
        m = length(clist)
    end
    if m != 2
        return false
    end

    # check that the two half-spaces are complementary
    return @inbounds iscomplement(clist[1], clist[2])
end

function load_symbolics_hpolyhedron()
    return quote
        """
            HPolyhedron(expr::Vector{<:Num}, vars=_get_variables(expr);
                        N::Type{<:Real}=Float64)

        Return a polyhedron in constraint representation given by a list of symbolic
        expressions.

        ### Input

        - `expr` -- vector of symbolic expressions that describes each half-space
        - `vars` -- (optional, default: `_get_variables(expr)`), if an array of
                    variables is given, use those as the ambient variables in the set
                    with respect to which derivations take place; otherwise, use only
                    the variables that appear in the given expression (but be careful
                    because the order may be incorrect; it is advised to always pass
                    `vars` explicitly)
        - `N`    -- (optional, default: `Float64`) the numeric type of the returned set

        ### Output

        An `HPolyhedron`.

        ### Examples

        ```jldoctest
        julia> using Symbolics

        julia> vars = @variables x y
        2-element Vector{Num}:
         x
         y

        julia> HPolyhedron([x + y <= 1, x + y >= -1], vars)
        HPolyhedron{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0), HalfSpace{Float64, Vector{Float64}}([-1.0, -1.0], 1.0)])

        julia> X = HPolyhedron([x == 0, y <= 0], vars)
        HPolyhedron{Float64, Vector{Float64}}(HalfSpace{Float64, Vector{Float64}}[HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], -0.0), HalfSpace{Float64, Vector{Float64}}([-1.0, -0.0], 0.0), HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], -0.0)])
        ```
        """
        function HPolyhedron(expr::Vector{<:Num}, vars::AbstractVector{Num};
                             N::Type{<:Real}=Float64)
            clist = Vector{HalfSpace{N,Vector{N}}}()
            sizehint!(clist, length(expr))
            got_hyperplane = false
            got_halfspace = false
            zeroed_vars = Dict(v => zero(N) for v in vars)
            vars_list = collect(vars)
            for ex in expr
                exval = Symbolics.value(ex)
                got_hyperplane, sexpr = _is_hyperplane(exval)
                if !got_hyperplane
                    got_halfspace, sexpr = _is_halfspace(exval)
                    if !got_halfspace
                        throw(ArgumentError("expected an expression describing either " *
                                            "a half-space of a hyperplane, got $expr"))
                    end
                end

                coeffs = [N(α.val) for α in gradient(sexpr, vars_list)]
                β = -N(Symbolics.substitute(sexpr, zeroed_vars))

                push!(clist, HalfSpace(coeffs, β))
                if got_hyperplane
                    push!(clist, HalfSpace(-coeffs, -β))
                end
            end
            return HPolyhedron(clist)
        end

        function HPolyhedron(expr::Vector{<:Num}; N::Type{<:Real}=Float64)
            return HPolyhedron(expr, _get_variables(expr); N=N)
        end
        function HPolyhedron(expr::Vector{<:Num}, vars; N::Type{<:Real}=Float64)
            return HPolyhedron(expr, _vec(vars); N=N)
        end
    end
end  # quote / load_symbolics_hpolyhedron()
