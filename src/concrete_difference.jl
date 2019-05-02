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
    difference(I1::Interval{N}, I2::Interval{N}) where {N}

Return the set difference between the given intervals.

The set difference is defined as:

```math
    I1 \\setminus I2 = \\{x: x ∈ I₁ \\text{ and } x ∉ I₂ \\}
```

### Input

- `I1` -- first interval
- `I2` -- second interval

### Output

- An `EmptySet` if the set difference is empty.
- A `UnionSetArray` of `Interval` sets otherwise.
"""
function difference(I1::Interval{N}, I2::Interval{N}) where {N}
    I12 = intersection(I1, I2)
    if isempty(I12)
        Idiff = [I1]
    else
        Ileft = Interval(min(I1), min(I12))
        Iright = Interval(max(I12), max(I1))

        Ileft_iszero = isapproxzero(IntervalArithmetic.diam(Ileft.dat))
        Iright_iszero = isapproxzero(IntervalArithmetic.diam(Iright.dat))

        if Ileft_iszero && Iright_iszero
            return EmptySet{N}()
        elseif Ileft_iszero && !Iright_iszero
            Idiff = [Iright]
        elseif !Ileft_iszero && Iright_iszero
            Idiff = [Ileft]
        else
            Idiff = [Ileft, Iright]
        end
    end
    return UnionSetArray(Idiff)
end
