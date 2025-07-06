"""
# Extended help

    rectify(X::Interval)

### Output

An `Interval`, even if it represents a singleton only containing the origin
(which is the case if the original interval was nonpositive).
"""
function rectify(X::Interval)
    N = eltype(X)
    if X.dat.lo >= zero(N)
        # interval is already nonnegative
        return X
    else
        # lower end is negative
        return Interval(zero(N), Base.max(X.dat.hi, zero(N)))
    end
end
