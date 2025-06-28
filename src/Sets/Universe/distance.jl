# distance point <-> set

@commutative function distance(x::AbstractVector, U::Universe; p::Real=2)
    @assert length(x) == dim(U) "incompatible dimensions $(length(x)) and $(dim(U))"

    N = promote_type(eltype(x), eltype(U))
    return N(0)
end

# distance set <-> set

function distance(U1::Universe, U2::Universe; p::Real=2.0)
    return _distance_universe(U1, U2; p=p)
end

function _distance_universe(U::Universe, X::LazySet; p::Real=2.0)
    @assert dim(U) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(U)) and $(dim(X)), respectively"

    N = promote_type(eltype(U), eltype(X))
    if isempty(X)
        return N(Inf)
    end
    return zero(N)
end
