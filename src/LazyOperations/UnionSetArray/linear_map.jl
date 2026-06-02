@validate function linear_map(M::AbstractMatrix, cup::UnionSetArray)
    return UnionSetArray([linear_map(M, X) for X in cup])
end
