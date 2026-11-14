function triangulate(::Universe; kwargs...)
    return throw(ArgumentError("cannot compute a triangulation for a universal set"))
end
