function project(X::Interval, block::AbstractVector{Int}; kwargs...)
    @assert length(block) == 1 && block[1] == 1 "invalid permutation vector $block"
    return X
end
