function scale(α::Real, Z::Zonotope)
    return _scale_copy_inplace(α, Z)
end

"""
# Extended help

    scale!(α::Real, Z::Zonotope)

### Algorithm

The result is obtained by applying the numerical scale to the center and
generators.
"""
function scale!(α::Real, Z::Zonotope)
    Z.center .*= α
    Z.generators .*= α
    return Z
end
