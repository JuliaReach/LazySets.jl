function σ(d::AbstractVector, X::Star)
    A = basis(X)
    return A * σ(At_mul_B(A, d), predicate(X)) + center(X)
end
