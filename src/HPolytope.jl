export HPolytope

"""
    HPolytope{N<:Real} <: AbstractPolytope{N}

Type that represents a convex polytope in H-representation.

### Fields

- `constraints` -- vector of linear constraints

### Note

Recall that a polytope is a bounded polyhedron. Boundedness is a running
assumption in this type.
"""
struct HPolytope{N<:Real} <: AbstractPolytope{N}
    constraints::Vector{LinearConstraint{N}}
end

# constructor for an HPolytope with no constraints
HPolytope{N}() where {N<:Real} = HPolytope{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolytope with no constraints of type Float64
HPolytope() = HPolytope{Float64}()

# conversion constructor
HPolytope(S::LazySet) = convert(HPolytope, S)

# constructor for an HPolytope from a simple H-representation
function HPolytope(A::Matrix{N}, b::Vector{N}) where {N<:Real}
    m = size(A, 1)
    constraints = LinearConstraint{N}[]
    @inbounds for i in 1:m
        push!(constraints, LinearConstraint(A[i, :], b[i]))
    end
    return HPolytope(constraints)
end

# ==========================================
# Lower level methods that use Polyhedra.jl
# ==========================================

function load_polyhedra_hpolytope() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

function convert(::Type{HPolytope{N}}, P::HRep{T, N}) where {T, N}
    constraints = LinearConstraint{N}[]
    for hi in Polyhedra.allhalfspaces(P)
        push!(constraints, HalfSpace(hi.a, hi.β))
    end
    return HPolytope(constraints)
end

"""
    HPolytope(P::HRep{T, N}, backend=default_polyhedra_backend(N)) where {T, N}

Return a polytope in H-representation given a `HRep` polyhedron
from `Polyhedra.jl`.

### Input

- `P` -- `HRep` polyhedron

### Output

An `HPolytope`.
"""
function HPolytope(P::HRep{T, N}) where {T, N}
    convert(HPolytope, P)
end

end # quote
end # function load_polyhedra_hpolytope()
