function permute(U::Universe, p::AbstractVector{Int})
    @assert length(p) == dim(U) "the dimensions should match, but they are $(length(p)) and " *
                                "$(dim(U)), respectively"
    @assert all(1 <= v <= dim(U) for v in p) "invalid dimension in index vector"

    return U
end
