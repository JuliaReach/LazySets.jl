"""
    difference(X::AbstractHyperrectangle, Y::AbstractHyperrectangle)

Compute the set difference between two hyperrectangular sets.

### Input

- `X` -- first hyperrectangular set
- `Y` -- second hyperrectangular set

The set difference is defined as:

```math
    X ∖ Y = \\{x: x ∈ X \\text{ and } x ∉ Y \\}
```

### Output

A `UnionSetArray` consisting of the union of hyperrectangles. Note that this
union is in general not convex.

### Algorithm

This implementation uses `IntervalArithmetic.setdiff`.
"""
function difference(X::AbstractHyperrectangle, Y::AbstractHyperrectangle)
    Xib = convert(IA.IntervalBox, X)
    Yib = convert(IA.IntervalBox, Y)
    return UnionSetArray(convert.(Hyperrectangle, IA.setdiff(Xib, Yib)))
end

function difference(X::Interval{N}, H::HalfSpace) where {N}
    @assert dim(H) == 1

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

function difference(X::LazySet{N}, U::Universe) where {N}
    n = dim(X)
    @assert n == dim(U) "incompatible dimensions $n and $(dim(U))"
    return EmptySet{N}(n)
end
