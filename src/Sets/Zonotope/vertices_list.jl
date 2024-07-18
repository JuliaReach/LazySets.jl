function _vertices_list_2D(c::AbstractVector{N}, G::AbstractMatrix{N};
                           apply_convex_hull::Bool) where {N}
    if apply_convex_hull
        return _vertices_list_zonotope_iterative(c, G; apply_convex_hull=apply_convex_hull)
    end

    require(@__MODULE__, :LazySets; fun_name="_vertices_list_2D")

    angles = mapslices(_angles, G; dims=1)[1, :]
    perm = sortperm(angles)
    sorted_angles = angles[perm]
    sorted_G = G[:, perm]
    polygons = Vector{VPolygon{N}}()
    sizehint!(polygons, 4)

    @inbounds for i in zip(0:90:360, 90:90:360)
        index = i[1] .<= sorted_angles .< i[2]
        if sum(index) > 0
            push!(polygons, _single_quadrant_vertices_enum(sorted_G[:, index], true))
        end
    end

    return vertices_list(translate(reduce(minkowski_sum, polygons), c))
end

@inline function _angles(point::AbstractVector{N}) where {N}
    return (atand(point[2], point[1]) + 360) % 360
end

function _single_quadrant_vertices_enum(G::AbstractMatrix{N}, sorted::Bool=false) where {N}
    if !sorted
        G = sortslices(G; dims=2, by=_angles)
    end
    return VPolygon(2 * cumsum(hcat(G, -G); dims=2) .- sum(G; dims=2))
end
