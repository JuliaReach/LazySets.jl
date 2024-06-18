"""
    rectify(x::Interval{N}) where {N}

Concrete rectification of an interval.

### Input

- `x` -- interval

### Output

The `Interval` that corresponds to the rectification of `x`.

### Notes

Note that the result is an `Interval` even if the set becomes a singleton (which
is the case if the original interval was nonpositive).
"""
function rectify(x::Interval{N}) where {N}
    if x.dat.lo >= zero(N)
        # interval is already nonnegative
        return x
    else
        # lower end is negative
        return Interval(zero(N), Base.max(x.dat.hi, zero(N)))
    end
end
