function constraints_list(L::Line2D)
    return _constraints_list_hyperplane(L.a, L.b)
end
