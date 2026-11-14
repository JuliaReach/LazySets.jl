function triangulate(::EmptySet; kwargs...)
    return throw(ArgumentError("cannot compute a triangulation for an empty set"))
end
