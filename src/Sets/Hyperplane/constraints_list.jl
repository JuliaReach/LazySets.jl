function constraints_list(H::Hyperplane)
    return _constraints_list_hyperplane(H.a, H.b)
end

# see ext/LazySets/LazySetsHyperplaneExt.jl
_constraints_list_hyperplane(a, b) = error()
