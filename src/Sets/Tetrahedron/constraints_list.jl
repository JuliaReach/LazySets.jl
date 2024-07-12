function constraints_list(T::Tetrahedron)
    return constraints_list(convert(VPolytope, T))
end
