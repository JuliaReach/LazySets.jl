function convex_hull(P::Polygon)
    require(@__MODULE__, :LazySets; fun_name="convex_hull")

    return VPolygon(P.vertices)
end
