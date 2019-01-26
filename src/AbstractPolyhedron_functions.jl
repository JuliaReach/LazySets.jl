import Base.∈

export constrained_dimensions

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
