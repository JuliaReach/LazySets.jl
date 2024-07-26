function Base.convert(::Type{VPolytope}, T::Tetrahedron)
    return VPolytope(T.vertices)
end
