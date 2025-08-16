@validate function low(Z::Zonotope, i::Int)
    G = genmat(Z)
    v = center(Z, i)
    @inbounds for j in 1:ngens(Z)
        v -= abs(G[i, j])
    end
    return v
end
