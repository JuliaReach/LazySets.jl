"""
    σ(d::AbstractVector, T::Tetrahedron)

Return a support vector of a tetrahedron in a given direction.

### Input

- `d` -- direction
- `T` -- tetrahedron

### Output

A support vector in the given direction.

### Algorithm

Currently falls back to the `VPolytope` implementation.
"""
function σ(d::AbstractVector, T::Tetrahedron)
    require(@__MODULE__, :LazySets; fun_name="σ")

    return σ(d, convert(VPolytope, T))
end
