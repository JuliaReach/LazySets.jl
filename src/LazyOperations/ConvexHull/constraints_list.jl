function constraints_list(ch::ConvexHull)
    ST = (dim(ch) == 2) ? VPolygon : VPolytope
    V = convert(ST, ch)
    return constraints_list(V)
end
