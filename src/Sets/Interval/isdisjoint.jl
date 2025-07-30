"""
# Extended help

    isdisjoint(I1::Interval, I2::Interval, [witness]::Bool=false)

### Algorithm

``I1 ∩ I2 = ∅`` iff there is a gap between the left-most point of the second
interval and the right-most point of the first interval, or vice-versa.

A witness is computed by taking the maximum over the left-most points of each
interval, which is guaranteed to belong to the intersection.
"""
@validate function isdisjoint(I1::Interval, I2::Interval, witness::Bool=false)
    if witness
        return _isdisjoint(I1, I2, Val(true))
    else
        return _isdisjoint(I1, I2, Val(false))
    end
end

function _isdisjoint(I1::Interval, I2::Interval, ::Val{false})
    return !_leq(min(I2), max(I1)) || !_leq(min(I1), max(I2))
end

function _isdisjoint(I1::Interval, I2::Interval, ::Val{true})
    if _isdisjoint(I1, I2, Val(false))
        N = promote_type(eltype(I1), eltype(I2))
        return (true, N[])
    else
        # mid-point, most robust witness
        l = max(min(I1), min(I2))
        h = min(max(I1), max(I2))
        w = [l + (h - l) / 2]
        return (false, w)
    end
end
