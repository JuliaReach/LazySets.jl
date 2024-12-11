function vertices(∅::EmptySet)
    N = eltype(∅)
    return EmptyIterator{Vector{N}}()
end
