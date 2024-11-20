"""
# Extended help

    rectify(x::Interval{N}) where {N}

### Output

An `Interval`, even if it represents a singleton only containing the origin
(which is the case if the original interval was nonpositive).
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
