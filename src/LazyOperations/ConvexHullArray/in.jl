# membership in a convex-hull array of singletons
@validate function in(x::AbstractVector, X::ConvexHullArray)
    n = length(x)
    ST = n == 2 ? VPolygon : VPolytope
    V = convert(ST, X)
    return x ∈ V
end
