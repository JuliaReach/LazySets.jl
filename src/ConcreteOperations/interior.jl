export interior


"""
    interior(d::AbstractVector, P::LazySet; ε=_TOL_F64.rtol)

Check if the point `d` is contained in the interior of the convex set `P`.

### Algorithm

The implementation checks if a `BallInf` with center `d` and radius `ε` is
contained in the set `P` which results in a numerical check for `d ∈ interior(P)`
with error tolerance `ε` with default value `_TOL_F64.rtol`.
"""
function interior(d::AbstractVector, P::LazySet; ε=_TOL_F64.rtol)
    return BallInf(d, ε) ⊆ P
end
