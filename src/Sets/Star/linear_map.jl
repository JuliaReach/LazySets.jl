function linear_map(M::AbstractMatrix, X::Star)
    c′ = M * X.c
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end
