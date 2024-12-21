"""
# Extended help

    ∈(v::AbstractVector, X::Star)

### Algorithm

See [`∈(::AbstractVector, ::LazySets.AbstractAffineMap)`](@ref).
"""
function ∈(x::AbstractVector, X::Star)
    return basis(X) \ (x - center(X)) ∈ predicate(X)
end
