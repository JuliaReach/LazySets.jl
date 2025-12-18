"""
# Extended help

    difference(X::AbstractHyperrectangle, Y::AbstractHyperrectangle)

### Output

A `UnionSetArray` consisting of a union of hyperrectangles describing the
difference, or an `EmptySet` if the difference is empty. Note that the union is
in general not convex.

### Algorithm

This implementation uses `IntervalArithmetic.setdiff`.
"""
@validate function difference(X::AbstractHyperrectangle, Y::AbstractHyperrectangle)
    Xib = convert(IA.IntervalBox, X)
    Yib = convert(IA.IntervalBox, Y)
    U = IA.setdiff(Xib, Yib)
    if isempty(U)
        return EmptySet(dim(X))
    end
    return UnionSetArray(convert.(Hyperrectangle, U))
end

@validate function difference(X::Interval{N}, H::HalfSpace) where {N}
    if H.a[1] < zero(N)
        # half-space is a lower bound
        l = low(H, 1)
        if l > max(X)
            # no effect
            return X
        elseif l <= min(X)
            # empty difference
            return EmptySet{N}(1)
        else
            return Interval(min(X), l)
        end
    else
        # half-space is an upper bound
        h = high(H, 1)
        if h < min(X)
            # no effect
            return X
        elseif h >= max(X)
            # empty difference
            return EmptySet{N}(1)
        else
            return Interval(h, max(X))
        end
    end
end

@validate function difference(U::Universe, X::LazySet)
    return _difference_universe(U, X)
end

@validate function difference(X::LazySet, U::Universe)
    return _difference_universe2(X, U)
end

@validate function difference(∅::EmptySet, X::LazySet)
    return _difference_emptyset(∅, X)
end

@validate function difference(X::LazySet, ∅::EmptySet)
    return _difference_emptyset2(X, ∅)
end

# ============== #
# disambiguation #
# ============== #

@validate function difference(∅::EmptySet, U::Universe)
    return _difference_emptyset(∅, U)
end

@validate function difference(U::Universe, ∅::EmptySet)
    return _difference_emptyset2(U, ∅)
end
