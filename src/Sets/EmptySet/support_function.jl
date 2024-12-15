function ρ(d::AbstractVector, ∅::EmptySet)
    @assert length(d) == dim(∅) "incompatible dimensions $(length(d)) and $(dim(∅))"

    throw(ArgumentError("the support function of an empty set is undefined"))
end
