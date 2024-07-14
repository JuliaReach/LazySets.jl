function permute(P::VPolytope, p::AbstractVector{Int})
    vlist = [v[p] for v in P.vertices]
    return VPolytope(vlist)
end
