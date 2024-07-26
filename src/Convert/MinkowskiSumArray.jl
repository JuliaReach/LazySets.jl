"""
    convert(::Type{MinkowskiSumArray},
            X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}

Convert the Minkowski sum of a Minkowski sum array to a Minkowski sum array.

### Input

- `MinkowskiSumArray`  -- target type
- `X`                  -- Minkowski sum of a Minkowski sum array

### Output

A Minkowski sum array.
"""
function convert(::Type{MinkowskiSumArray},
                 X::MinkowskiSum{N,ST,MinkowskiSumArray{N,ST}}) where {N,ST}
    return MinkowskiSumArray(vcat(first(X), array(second(X))))
end
