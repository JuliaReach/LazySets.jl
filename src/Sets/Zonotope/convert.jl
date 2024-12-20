function convert(::Type{Zonotope}, Z::AbstractZonotope)
    return Zonotope(center(Z), genmat(Z))
end
