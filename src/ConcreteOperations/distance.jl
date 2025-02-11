function distance(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle;
                  p::Real=2)
    n = dim(H1)
    @assert n == dim(H2) "incompatible set dimensions $n and $(dim(H2))"

    N = promote_type(eltype(H1), eltype(H2))
    d = Vector{N}(undef, n)
    @inbounds for i in 1:n
        # find relative position in dimension i
        # (if c1 == c2, the results are equivalent independent of the branch)
        if center(H1, i) >= center(H2, i)
            lhs = low(H1, i)
            rhs = high(H2, i)
        else
            lhs = low(H2, i)
            rhs = high(H1, i)
        end
        if _leq(lhs, rhs)
            d[i] = zero(N)
        else
            d[i] = rhs - lhs
        end
    end
    return norm(d, p)
end

@commutative function distance(S::AbstractSingleton, X::LazySet; p::Real=2.0)
    return distance(element(S), X; p=p)
end

@commutative function distance(∅::EmptySet, X::LazySet; p::Real=2.0)
    return _distance_emptyset(∅, X; p=p)
end

@commutative function distance(U::Universe, X::LazySet; p::Real=2.0)
    return _distance_universe(U, X; p=p)
end

# ============== #
# disambiguation #
# ============== #

function distance(S1::AbstractSingleton, S2::AbstractSingleton; p::Real=2.0)
    return distance(element(S1), element(S2); p=p)
end

@commutative function distance(S::AbstractSingleton, H::AbstractHyperrectangle; p::Real=2.0)
    return distance(element(S), H; p=p)
end

@commutative function distance(U::Universe, S::AbstractSingleton; p::Real=2.0)
    return _distance_universe(U, S; p=p)
end

for T in [:AbstractHyperrectangle, :AbstractSingleton, :Universe]
    @eval @commutative function distance(∅::EmptySet, X::($T); p::Real=2.0)
        return _distance_emptyset(∅, X; p=p)
    end
end
