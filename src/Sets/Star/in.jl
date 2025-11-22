"""
# Extended help

    in(v::AbstractVector, X::Star)

### Algorithm

See [`in(::AbstractVector, ::LazySets.AbstractAffineMap)`](@ref).
"""
@validate function in(x::AbstractVector, X::Star)
    return basis(X) \ (x - center(X)) ∈ predicate(X)
end
