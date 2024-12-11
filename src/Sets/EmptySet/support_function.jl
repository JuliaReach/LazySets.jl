function ρ(::AbstractVector, ::EmptySet)
    throw(ArgumentError("the support function of an empty set is undefined"))
end
