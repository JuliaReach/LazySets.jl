function triangulate(::EmptySet; kwargs...)
    throw(ArgumentError("cannot compute a triangulation for an empty set"))
end
