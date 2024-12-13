function ∈(v::AbstractVector, ∅::EmptySet)
    @assert length(v) == dim(∅) "the dimensions should match, but they are $(length(v)) and " *
                                "$(dim(∅)), respectively"
    return false
end
