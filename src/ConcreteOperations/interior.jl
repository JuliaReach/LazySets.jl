export is_interior_point


"""
    is_interior_point(d::AbstractVector, P::LazySet; p=Inf, ε=_rtol(N))

Check if the point `d` is contained in the interior of the convex set `P`.

### Input

- `d`  -- point
- `P`  -- set
- `p`  -- norm of the ball used to apply the error tolerance
- `ε`  -- error tolerance of check

### Output

Boolean which indicates if point `d` is contained in `P`.

### Algorithm

The implementation checks if a `Ballp` of norm `p` with center `d` and radius `ε` is
contained in the set `P`. This is a numerical check for `d ∈ interior(P)`
with error tolerance `ε`.
"""
function is_interior_point(d::AbstractVector{N}, P::LazySet{N}; p=Inf, ε=_rtol(N)) where {N<:Real}
    return Ballp(p, d, ε) ⊆ P
end
