"""
    in(x::AbstractVector, R::Rectification)

Check whether a given point is contained in a rectification.

### Input

- `x` -- point/vector
- `R` -- rectification

### Output

`true` iff ``x ∈ R``.

### Algorithm

We first scan for negative entries in the vector.
If there are any, the vector is not contained in the rectification.

Next we ask a membership query in the wrapped set.
If the answer is positive, the vector is contained in the rectification.
(This holds because negative entries have been ruled out before.)

Otherwise, we scan for zero entries in the vector.
If there are none, membership reduces to membership in the wrapped set, and so
the answer is negative.

Finally, if there are zero entries in the vector and the vector is not contained
in the wrapped set, we give up and throw an error.
"""
@validate function in(x::AbstractVector, R::Rectification)
    N = promote_type(eltype(x), eltype(R))
    # scan for negative entries
    if any(xi -> xi < zero(N), x)
        return false
    end

    # membership test in the wrapped set
    if x ∈ R.X
        return true
    end

    # scan for zero entries
    if all(!iszero, x)
        return false
    end

    # compute exact set
    Y = _compute_exact_representation!(R)
    return x ∈ Y
end
