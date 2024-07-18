function permute(Z::Zonotope, p::AbstractVector{Int})
    c = Z.center[p]
    G = Z.generators[p, :]
    return Zonotope(c, G)
end
