function project(V::VPolygon, block::AbstractVector{Int}; kwargs...)
    if length(block) == 1
        @assert block[1] == 1 || block[1] == 2 "invalid projection to $block"
        l, h = extrema(V, block[1])
        return Interval(l, h)
    elseif length(block) == 2
        if block[1] == 1 && block[2] == 2
            return V  # no projection
        else
            @assert block[1] == 2 && block[2] == 1 "invalid projection to $block"
            return permute(V, block)  # swap dimensions
        end
    end
    throw(ArgumentError("invalid projection to $block"))
end
