@validate function constraints(U::Universe{N}) where {N}
    return EmptyIterator{Vector{N}}()
end
