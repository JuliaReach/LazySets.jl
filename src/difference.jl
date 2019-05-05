export difference

# alias for set difference
import Base: \

"""
    \\(X::LazySet, Y::LazySet)

Convenience alias for set difference.

### Input

- `X` -- a set
- `Y` -- another set

### Output

The set difference between `X` and `Y`.
"""
\(X::LazySet, Y::LazySet) = difference(X, Y)

# =================================
# Set difference between intervals
# =================================

"""
    difference(I1::IN, I2::IN)::Union{EmptySet{N}, IN, UnionSet{N, IN, IN}} where {N, IN<:Interval{N}}

Return the set difference between the given intervals.

The set difference is defined as:

```math
    I₁ \\setminus I₂ = \\{x: x ∈ I₁ \\text{ and } x ∉ I₂ \\}
```

### Input

- `I1` -- first interval
- `I2` -- second interval

### Output

Depending on the position of the intervals, the output is one of the following:

- An `EmptySet`.
- An `Interval`.
- A `UnionSet` of two `Interval` sets.

### Algorithm

Let ``I₁ = [a, b]`` and ``I₂ = [c, d]`` be intervals. Their set difference is
``I₁ \\setminus I₂ = \\{x: x ∈ I₁ \\text{ and } x ∉ I₂ \\}`` and depending on their
position three different results may occur:

- If ``I₁`` and ``I₂`` do not overlap, i.e. if their intersection is empty, then the
  set difference is just ``I₁``.
- Otherwise, let `I₁₂ = I₁ ∩ I₂` and assume that it is not empty, then either
  ``I₁₂`` splits `I₁` into one interval or into two intervals. The latter case
   happens when the inclusion is strict.

To check for strict inclusion, we asssume that the inclusion is strict and then
check if the resulting intervals that cover `I₁` (one to its left and one to its
right, let them be `Ileft` and `Iright`), obtained by intersection with `I₂`,
are flat or not. Three cases may arise:

- If both `Ileft` and `Iright` are flat then it means that `I₁ = I₂`, then the
  set difference is the empty set.
- If only `Ileft` is flat, then the remaining interval not covered by `I₂` is
  `Iright`. In a similar manner, if only `Iright` is flat, then `Ileft` is returned.
- Finally, if none of the intervals is flat, then `I₂` is strictly contained in
  `I₁` and the set union of `Ileft` and `Iright` is returned.
"""
function difference(I1::IN, I2::IN)::Union{EmptySet{N}, IN, UnionSet{N, IN, IN}} where {N, IN<:Interval{N}}
    I12 = intersection(I1, I2)
    if isempty(I12)
        return I1
    else
        Ileft = Interval(min(I1), min(I12))
        Iright = Interval(max(I12), max(I1))

        flat_left = isflat(Ileft)
        flat_right = isflat(Iright)

        if flat_left && flat_right
            return EmptySet{N}()
        elseif flat_left && !flat_right
            return Iright
        elseif !flat_left && flat_right
            return Ileft
        else
            return UnionSet(Ileft, Iright)
        end
    end
end

# return whether is flat, i.e. if the extreme values or the interval coincide.
function isflat(I::Interval)
    return isapproxzero(IntervalArithmetic.diam(I.dat))
end
