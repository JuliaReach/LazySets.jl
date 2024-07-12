"""
    σ(d::AbstractVector, P::VPolytope)

Return a support vector of a tetrahedron in a given direction.

### Input

- `d` -- direction
- `P` -- tetrahedron

### Output

A support vector in the given direction.

### Algorithm

Currently falls back to the `VPolytope` implementation.
"""
function σ(d::AbstractVector, T::Tetrahedron)
    return σ(d, convert(VPolytope, T))
end
