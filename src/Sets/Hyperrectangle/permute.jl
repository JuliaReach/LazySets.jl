function permute(H::Hyperrectangle, p::AbstractVector{Int})
    c = H.center[p]
    r = H.radius[p]
    return Hyperrectangle(c, r)
end
