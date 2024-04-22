export difference

# alias for set difference
import Base: \

"""
    \\(X::LazySet, Y::LazySet)

Convenience alias for set difference.

### Input

- `X` -- first set
- `Y` -- second set

### Output

The set difference between `X` and `Y`.

### Notes

If `X` and `Y` are intervals, `X \\ Y` is used in some libraries to denote
the left division, as the example below shows. However, it should not be
confused with the *set difference*. For example,

```jldoctest
julia> X = Interval(0, 2); Y = Interval(1, 4);

julia> X \\ Y  # computing the set difference
Interval{Float64}([0, 1])

julia> X.dat \\ Y.dat  # computing the left division
[0.5, ∞]
```
"""
\(X::LazySet, Y::LazySet) = difference(X, Y)

"""
    difference(X::Interval{N}, Y::Interval) where {N}

Compute the set difference between two intervals.

The set difference is defined as:

```math
    X \\setminus Y = \\{x: x ∈ X \\text{ and } x ∉ Y \\}
```

The backslash symbol, `\\`, can be used as an alias.

### Input

- `X` -- first interval
- `Y` -- second interval

### Output

Depending on the position of the intervals, the output is one of the following:

- An `EmptySet`.
- An `Interval`.
- A `UnionSet` of two `Interval` sets.

### Algorithm

Let ``X = [a, b]`` and ``Y = [c, d]`` be intervals. Their set difference is
``X ∖ Y = \\{x: x ∈ X \\text{ and } x ∉ Y \\}`` and, depending on their
position, three different results may occur:

- If ``X`` and ``Y`` do not overlap, i.e., if their intersection is empty, then
  the set difference is just ``X``.
- Otherwise, let ``Z = X ∩ Y ≠ ∅``, then ``Z`` splits ``X`` into either one or
  two intervals. The latter case happens when the bounds of ``Y`` are strictly
  contained in ``X``.

To check for strict inclusion, we assume that the inclusion is strict and then
check if the resulting intervals that cover ``X`` (one to its left and one to
its right, let them be `L` and `R`), obtained by intersection with ``Y``, are
flat or not. Three cases may arise:

- If both `L` and `R` are flat then ``X = Y`` and the result is the empty set.
- If only `L` is flat, then the result is `R`, the remaining interval not
  covered by ``Y``. Similarly, if only `R` is flat, then the result is `L`.
- Finally, if none of the intervals is flat, then ``Y`` is strictly contained
  in ``X`` and the set union of `L` and `R` is returned.
"""
function difference(X::Interval{N}, Y::Interval) where {N}
    Z = intersection(X, Y)
    if isempty(Z)
        return X
    else
        L = Interval(min(X), min(Z))
        R = Interval(max(Z), max(X))

        flat_left = isflat(L)
        flat_right = isflat(R)

        if flat_left && flat_right
            return EmptySet{N}(1)
        elseif flat_left && !flat_right
            return R
        elseif !flat_left && flat_right
            return L
        else
            return UnionSet(L, R)
        end
    end
end

"""
    difference(X::AbstractHyperrectangle{N}, Y::AbstractHyperrectangle) where {N}

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
function difference(X::AbstractHyperrectangle{N}, Y::AbstractHyperrectangle) where {N}
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
