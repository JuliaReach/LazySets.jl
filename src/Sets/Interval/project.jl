function project(X::Interval, block::AbstractVector{Int}; kwargs...)
    @assert length(block) == 1 && (@inbounds block[1]) == 1 "invalid permutation vector $block " *
                                                            "for Interval"
    return X
end
