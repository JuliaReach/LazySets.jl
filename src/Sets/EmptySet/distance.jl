# distance point <-> set

@validate_commutative function distance(x::AbstractVector, ∅::EmptySet; p::Real=2)
    N = promote_type(eltype(x), eltype(∅))
    return N(Inf)
end

# distance set <-> set

@validate function distance(∅₁::EmptySet, ∅₂::EmptySet; p::Real=2.0)
    return _distance_emptyset(∅₁, ∅₂; p=p)
end

function _distance_emptyset(∅::EmptySet, X::LazySet; p::Real=2.0)
    N = promote_type(eltype(∅), eltype(X))
    return N(Inf)
end
