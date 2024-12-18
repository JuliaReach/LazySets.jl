"""
# Extended help

    σ(d::AbstractVector, T::Tetrahedron)

### Algorithm

This method falls back to the `VPolytope` implementation.
"""
function σ(d::AbstractVector, T::Tetrahedron)
    require(@__MODULE__, :LazySets; fun_name="σ")

    return σ(d, convert(VPolytope, T))
end
