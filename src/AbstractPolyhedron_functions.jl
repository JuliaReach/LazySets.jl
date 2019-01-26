import Base.∈

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
