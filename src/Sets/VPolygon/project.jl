@validate function project(V::VPolygon, block::AbstractVector{Int}; kwargs...)
    if length(block) == 1
        l, h = extrema(V, block[1])
        return Interval(l, h)
    end
    # length(block) == 2
    if block[1] == 1 && block[2] == 2
        return V  # no projection
    else
        # block[1] == 2 && block[2] == 1
        return permute(V, block)  # swap dimensions
    end
end
