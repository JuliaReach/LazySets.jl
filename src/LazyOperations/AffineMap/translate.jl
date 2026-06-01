@validate function translate(am::AffineMap, x::AbstractVector)
    M = matrix(am)
    X = set(am)
    v = vector(am)
    return AffineMap(M, X, v + x)
end
