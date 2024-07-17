function permute(P::HPoly, p::AbstractVector{Int})
    T = basetype(P)
    return T([permute(H, p) for H in P.constraints])
end
