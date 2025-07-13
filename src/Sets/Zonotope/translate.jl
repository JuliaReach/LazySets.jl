@validate function translate!(Z::Zonotope, v::AbstractVector)
    Z.center .+= v
    return Z
end
