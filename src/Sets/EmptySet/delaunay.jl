function delaunay(::EmptySet)
    throw(ArgumentError("cannot compute a Delaunay triangulation for an empty set"))
end
