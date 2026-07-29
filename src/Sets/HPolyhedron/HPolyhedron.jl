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
