function constraints_list(T::Tetrahedron)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    return constraints_list(convert(VPolytope, T))
end
