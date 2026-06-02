function constraints_list(cha::ConvexHullArray)
    ST = (dim(cha) == 2) ? VPolygon : VPolytope
    V = convert(ST, cha)
    return constraints_list(V)
end
