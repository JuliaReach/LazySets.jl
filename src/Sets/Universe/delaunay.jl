function delaunay(::Universe)
    throw(ArgumentError("cannot compute a Delaunay triangulation for a universal set"))
end
