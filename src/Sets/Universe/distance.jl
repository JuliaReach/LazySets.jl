# distance point <-> set

@validate_commutative function distance(x::AbstractVector, U::Universe; p::Real=2)
    N = promote_type(eltype(x), eltype(U))
    return N(0)
end

# distance set <-> set

@validate function distance(U1::Universe, U2::Universe; p::Real=2.0)
    return _distance_universe(U1, U2; p=p)
end

function _distance_universe(U::Universe, X::LazySet; p::Real=2.0)
    N = promote_type(eltype(U), eltype(X))
    if isempty(X)
        return N(Inf)
    end
    return zero(N)
end
