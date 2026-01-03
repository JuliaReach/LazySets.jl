"""
    ballinf_approximation(S::LazySet)

Overapproximate a set by a tight ball in the infinity norm.

### Input

- `S` -- set

### Output

A tight ball in the infinity norm.

### Algorithm

The center and radius of the ball are obtained by averaging the low and high
coordinates of `S` computed with the `extrema` function.
"""
function ballinf_approximation(S::LazySet{N}) where {N}
    n = dim(S)
    c = Vector{N}(undef, n)
    r = zero(N)
    @inbounds for i in 1:n
        lo, hi = extrema(S, i)
        rcur = (hi - lo) / 2
        if (rcur > r)
            r = rcur
        elseif !_geq(rcur, zero(N))
            # contradicting bounds => set is empty
            throw(ArgumentError("set must be nonempty"))
        end
        c[i] = (hi + lo) / 2
    end
    return BallInf(c, r)
end
