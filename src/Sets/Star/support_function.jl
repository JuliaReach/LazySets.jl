function ρ(d::AbstractVector, X::Star)
    return ρ(At_mul_B(basis(X), d), predicate(X)) + dot(d, center(X))
end
