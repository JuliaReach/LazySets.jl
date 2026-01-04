@validate function translate(Z::Zonotope, v::AbstractVector)
    c = Z.center .+ v
    return Zonotope(c, Z.generators)
end

@validate function translate!(Z::Zonotope, v::AbstractVector)
    Z.center .+= v
    return Z
end
