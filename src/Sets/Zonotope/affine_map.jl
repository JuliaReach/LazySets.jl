@validate function affine_map(M, Z::Zonotope, v::AbstractVector)
    return translate!(linear_map(M, Z), v)
end
