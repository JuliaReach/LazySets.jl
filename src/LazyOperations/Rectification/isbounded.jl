"""
    isbounded(R::Rectification)

Check whether a rectification is bounded.

### Input

- `R` -- rectification

### Output

`true` iff the rectification is bounded.

### Algorithm

Let ``X`` be the set wrapped by rectification ``R``.
We first check whether ``X`` is bounded (because then ``R`` is bounded).
Otherwise, we check unboundedness of ``X`` in direction ``(1, 1, …, 1)``, which
is sufficient for unboundedness of ``R``; this step is not necessary but rather
a heuristics.
Otherwise, we check boundedness of ``X`` in every positive unit direction, which
is sufficient and necessary for boundedness of ``R``.
"""
function isbounded(R::Rectification)
    # check boundedness of the wrapped set
    if isbounded(R.X)
        return true
    end

    # heuristics: check boundedness in the one-vector direction
    n = dim(R.X)
    N = eltype(R)
    if isinf(ρ(ones(N, n), R.X))
        return false
    end

    # check boundedness in positive unit directions
    @inbounds for i in 1:n
        if isinf(high(R.X, i))
            return false
        end
    end
    return true
end
