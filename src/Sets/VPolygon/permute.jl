function permute(V::VPolygon, p::AbstractVector{Int})
    return VPolygon([v[p] for v in V.vertices]; apply_convex_hull=true)
end
