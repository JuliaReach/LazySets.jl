@validate function affine_map(M::AbstractMatrix, X::Star, v::AbstractVector)
    c′ = M * X.c + v
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end
