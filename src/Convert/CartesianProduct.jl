"""
    convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N}

Convert a two-dimensional hyperrectangle to the Cartesian product of two
intervals.

### Input

- `CartesianProduct` -- target type
- `H`                -- hyperrectangle

### Output

The Cartesian product of two intervals.
"""
function convert(::Type{CartesianProduct{N,Interval{N},Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N}
    @assert dim(H) == 2 "the hyperrectangle must be two-dimensional to " *
                        "convert it to the Cartesian product of two intervals, but it is " *
                        "$(dim(H))-dimensional; consider converting it to a " *
                        "`CartesianProductArray{$N, Interval{$N}}` instead"
    I1 = Interval(low(H, 1), high(H, 1))
    I2 = Interval(low(H, 2), high(H, 2))
    return CartesianProduct(I1, I2)
end
