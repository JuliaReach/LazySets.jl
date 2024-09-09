"""
    convert(::Type{CartesianProductArray{N, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N}

Convert a hyperrectangle to the Cartesian product array of intervals.

### Input

- `CartesianProductArray` -- target type
- `H`                     -- hyperrectangle

### Output

The Cartesian product of a finite number of intervals.
"""
function convert(::Type{CartesianProductArray{N,Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N}
    Iarray = [Interval(low(H, i), high(H, i)) for i in 1:dim(H)]
    return CartesianProductArray(Iarray)
end
