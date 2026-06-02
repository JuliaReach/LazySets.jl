function linear_map(M::AbstractMatrix, cup::UnionSet)
    return UnionSet(linear_map(M, cup.X), linear_map(M, cup.Y))
end
