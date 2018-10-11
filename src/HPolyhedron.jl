export HPolyhedron

"""
    HPolyhedron{N<:Real} <: LazySet{N}

Type that represents a convex polyhedron in H-representation.

### Fields

- `constraints` -- vector of linear constraints
"""
struct HPolyhedron{N<:Real} <: LazySet{N}
    constraints::Vector{LinearConstraint{N}}
end

# constructor for an HPolyhedron with no constraints
HPolyhedron{N}() where {N<:Real} = HPolyhedron{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolyhedron with no constraints of type Float64
HPolyhedron() = HPolyhedron{Float64}()

"""
    tosimplehrep(P::Union{HPolytope{N}, HPolyhedron{N}}) where {N}

Return the simple H-representation ``Ax ≤ b`` of a polyhedron.

### Input

- `P` -- polyhedron or polytope

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` are the offsets.
"""
function tosimplehrep(P::Union{HPolytope{N}, HPolyhedron{N}}) where {N}
    if length(P.constraints) == 0
        A = Matrix{N}(undef, 0, 0)
        b = Vector{N}(undef, 0)
        return (A, b)
    end

    A = hcat([ci.a for ci in P.constraints]...)'
    b = [ci.b for ci in P.constraints]
    return (A, b)
end