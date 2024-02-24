using .TaylorModels: domain

# check that a vector of Taylor models has the [-1, 1] domain
function _has_normalized_domain(vTM)
    return all(p -> all(==(IA.interval(-1, 1)), domain(p)), vTM)
end

eval(load_taylormodels_convert_polynomial_zonotope())
