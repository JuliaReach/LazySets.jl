function project(P::VPolytope, block::AbstractVector{Int}; kwargs...)
    require(@__MODULE__, :LazySets; fun_name="project")

    if isempty(P.vertices)
        return P
    end

    m = length(block)
    if m == 1
        l, h = extrema(P, block[1])
        return Interval(l, h)
    end

    M = projection_matrix(block, dim(P), eltype(P.vertices))
    πvertices = broadcast(v -> M * v, P.vertices)

    if m == 2
        return VPolygon(πvertices; apply_convex_hull=true)
    else
        return VPolytope(πvertices)
    end
end
