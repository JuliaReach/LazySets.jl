function translate(∅::EmptySet, v::AbstractVector)
    return translate!(∅, v)  # no need to copy
end

function translate!(∅::EmptySet, v::AbstractVector)
    @assert length(v) == dim(∅) "cannot translate a $(dim(∅))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return ∅
end
