@validate function σ(d::AbstractVector, ∅::EmptySet)
    throw(ArgumentError("the support vector of an empty set is undefined"))
end
