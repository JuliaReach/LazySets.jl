"""
    project(S::LazySet, block::AbstractVector{Int}, set_type::Type{<:LinearMap},
            [n]::Int=dim(S); [kwargs...])

Project a high-dimensional set to a given block by using a lazy linear map.

### Input

- `S`         -- set
- `block`     -- block structure - a vector with the dimensions of interest
- `LinearMap` -- used for dispatch
- `n`         -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A lazy `LinearMap` representing the projection of the set `S` to block `block`.
"""
@inline function project(S::LazySet, block::AbstractVector{Int},
                         set_type::Type{<:LinearMap}, n::Int=dim(S);
                         kwargs...)
    N = eltype(S)
    M = projection_matrix(block, n, N)
    return M * S
end
