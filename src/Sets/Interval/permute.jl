function permute(X::Interval, p::AbstractVector{Int})
    @assert length(p) == 1 && p[1] == 1 "invalid permutation vector $p"
    return X
end
