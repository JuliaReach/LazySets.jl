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
    I1 \\setminus I2 = \\{x: x ∈ I₁ \\text{ and } x ∉ I₂ \\}
```

### Input

- `I1` -- first interval
- `I2` -- second interval

### Output

Depending on the position of the intervals, the output is one of the following:

- An `EmptySet`.
- An `Interval`.
- A `UnionSet` of two `Interval` sets.
"""
function difference(I1::IN, I2::IN)::Union{EmptySet{N}, IN, UnionSet{N, IN, IN}} where {N, IN<:Interval{N}}
    I12 = intersection(I1, I2)
    if isempty(I12)
        return I1
    else
        Ileft = Interval(min(I1), min(I12))
        Iright = Interval(max(I12), max(I1))

        Ileft_iszero = isapproxzero(IntervalArithmetic.diam(Ileft.dat))
        Iright_iszero = isapproxzero(IntervalArithmetic.diam(Iright.dat))

        if Ileft_iszero && Iright_iszero
            return EmptySet{N}()
        elseif Ileft_iszero && !Iright_iszero
            return Iright
        elseif !Ileft_iszero && Iright_iszero
            return Ileft
        else
            return UnionSet(Ileft, Iright)
        end
    end
end
