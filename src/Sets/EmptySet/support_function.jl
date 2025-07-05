@validate function ρ(d::AbstractVector, ∅::EmptySet)
    throw(ArgumentError("the support function of an empty set is undefined"))
end
