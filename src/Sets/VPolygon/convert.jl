function convert(::Type{VPolygon}, X::LazySet)
    @assert dim(X) == 2 "set must be two-dimensional for conversion"
    return VPolygon(vertices_list(X))
end
