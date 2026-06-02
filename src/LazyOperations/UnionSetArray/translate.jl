@validate function translate(cup::UnionSetArray, v::AbstractVector)
    return UnionSetArray([translate(X, v) for X in cup])
end
