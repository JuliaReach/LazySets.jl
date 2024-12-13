function permute(∅::EmptySet, p::AbstractVector{Int})
    @assert length(p) == dim(∅) "the dimensions should match, but they are $(length(p)) and " *
                                "$(dim(∅)), respectively"
    @assert all(1 <= v <= dim(∅) for v in p) "invalid dimension in index vector"

    return ∅
end
