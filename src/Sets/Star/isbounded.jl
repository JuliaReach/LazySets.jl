"""
# Extended help

    isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)

### Algorithm

See [`isbounded(::LazySets.AbstractAffineMap)`](@ref).
"""
function isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)
    am = convert(STAR, X)
    return isbounded(am; cond_tol=cond_tol)
end
