function copy(Z::Zonotope)
    return Zonotope(copy(Z.center), copy(Z.generators))
end
