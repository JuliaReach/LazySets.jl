function permute(X::Interval, p::AbstractVector{Int})
    @assert length(p) == 1 && (@inbounds p[1]) == 1 "invalid permutation vector $p for Interval"
    return X
end
