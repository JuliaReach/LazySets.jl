@validate function translate(∅::EmptySet, v::AbstractVector)
    return translate!(∅, v)  # no need to copy
end

@validate function translate!(∅::EmptySet, v::AbstractVector)
    return ∅
end
