"""
    is_interior_point(v::AbstractVector{<:Real}, X::LazySet; kwargs...) end

Check whether a point is contained in the interior of a set.

### Input

- `v`  -- point/vector
- `X`  -- set
- `p`  -- (optional; default: `Inf`) norm of the ball used to apply the error
          tolerance
- `ε`  -- (optional; default: `_rtol(eltype(X))`) error tolerance of the check

### Output

`true` iff the point `v` is strictly contained in `X` with tolerance `ε`.
"""
function is_interior_point(::AbstractVector{<:Real}, ::LazySet; kwargs...) end
