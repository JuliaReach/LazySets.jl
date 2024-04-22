"""
    is_interior_point(v::AbstractVector{N}, X::LazySet; [p]=N(Inf), [ε]=_rtol(N)) where {N}

Check whether a point is contained in the interior of a set.

### Input

- `v`  -- point/vector
- `X`  -- set
- `p`  -- (optional; default: `N(Inf)`) norm of the ball used to apply the error
          tolerance
- `ε`  -- (optional; default: `_rtol(N)`) error tolerance of the check

### Output

`true` iff the point `v` is strictly contained in `X` with tolerance `ε`.
"""
function is_interior_point(::AbstractVector{N}, ::LazySet; p=N(Inf), ε=_rtol(N)) where {N} end
