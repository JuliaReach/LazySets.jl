function convex_hull(P::VPolygon, Q::VPolygon;
                     algorithm::String="monotone_chain")
    vunion = [P.vertices; Q.vertices]
    convex_hull!(vunion; algorithm=algorithm)
    return VPolygon(vunion; apply_convex_hull=false)
end
