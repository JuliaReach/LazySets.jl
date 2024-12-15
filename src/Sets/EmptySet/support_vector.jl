function σ(d::AbstractVector, ∅::EmptySet)
    @assert length(d) == dim(∅) "incompatible dimensions $(length(d)) and $(dim(∅))"

    throw(ArgumentError("the support vector of an empty set is undefined"))
end
