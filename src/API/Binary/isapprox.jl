"""
    ≈(X::LazySet, Y::LazySet)

Check whether two sets of the same type are approximately equal.

### Input

- `X` -- set
- `Y` -- set of the same base type as `X`

### Output

`true` iff `X` is approximately equal to `Y`.

### Notes

The check is purely syntactic and the sets need to have the same base type, i.e.,
`X::T1 ≈ Y::T2` always returns `false` even if `X` and `Y` represent the same
set. But `X::T{Int64} ≈ Y::T{Float64}` is a valid comparison. Note that, unlike
most other binary operations, a query with `X` and `Y` of different dimension is
allowed (and results in the answer `false`).

"`≈`" can be typed by `\\approx<tab>`.
The convenience alias `isapprox` is also available.
"""
function ≈(::LazySet, ::LazySet) end
