import Base.∈

export constrained_dimensions,
       tosimplehrep

"""
    ∈(x::AbstractVector{N}, P::AbstractPolyhedron{N})::Bool where {N<:Real}

Check whether a given point is contained in a polyhedron.

### Input

- `x` -- point/vector
- `P` -- polyhedron

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies inside each defining half-space.
"""
function ∈(x::AbstractVector{N}, P::AbstractPolyhedron{N})::Bool where {N<:Real}
    @assert length(x) == dim(P) "a $(length(x))-dimensional point cannot be " *
        "an element of a $(dim(P))-dimensional set"

    for c in constraints_list(P)
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end

"""
    constrained_dimensions(P::AbstractPolyhedron{N})::Vector{Int}
        where {N<:Real}

Return the indices in which a polyhedron is constrained.

### Input

- `P` -- polyhedron

### Output

A vector of ascending indices `i` such that the polyhedron is constrained in
dimension `i`.

### Examples

A 2D polyhedron with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(P::AbstractPolyhedron{N}
                               )::Vector{Int} where {N<:Real}
    zero_indices = zeros(Int, dim(P))
    for constraint in constraints_list(P)
        for i in constrained_dimensions(constraint)
            zero_indices[i] = i
        end
    end
    return filter(x -> x != 0, zero_indices)
end

"""
    tosimplehrep(constraints::AbstractVector{LinearConstraint{N}})
        where {N<:Real}

Return the simple H-representation ``Ax ≤ b`` from a list of linear constraints.

### Input

- `constraints` -- a list of linear constraints

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.
"""
function tosimplehrep(constraints::AbstractVector{LinearConstraint{N}}
                     ) where {N<:Real}
    n = length(constraints)
    if n == 0
        A = Matrix{N}(undef, 0, 0)
        b = Vector{N}(undef, 0)
        return (A, b)
    end
    A = zeros(N, n, dim(first(constraints)))
    b = zeros(N, n)
    @inbounds begin
        for (i, Pi) in enumerate(constraints)
            A[i, :] = Pi.a
            b[i] = Pi.b
        end
    end
    return (A, b)
end
