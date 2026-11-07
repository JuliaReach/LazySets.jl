@doc """
    project(X::LazySet, block::AbstractVector{Int})

Project a set to a given block by using a concrete linear map.

### Input

- `X`       -- set
- `block`   -- block structure - a vector with the dimensions of interest

### Output

A set representing the projection of `X` to block `block`.
"""
project(::LazySet, ::AbstractVector{Int})
function project end
