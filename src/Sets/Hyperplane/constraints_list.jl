@validate function constraints_list(H::Hyperplane)
    return _constraints_list_hyperplane(H.a, H.b)
end

# internal helper function
function _constraints_list_hyperplane(a::AbstractVector, b)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    return [HalfSpace(a, b), HalfSpace(-a, -b)]
end
