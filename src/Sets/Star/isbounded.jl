"""
    isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)

Check whether a star is bounded.

### Input

- `X`        -- star
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the basis matrix is invertible)

### Output

`true` iff the star is bounded.

### Algorithm

See [`isbounded(::LazySets.AbstractAffineMap)`](@ref).
"""
function isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)
    am = convert(STAR, X)
    return isbounded(am; cond_tol=cond_tol)
end
