"""
    isdisjoint(I1::Interval, I2::Interval, [witness]::Bool=false)

Check whether two intervals do not intersect, and otherwise optionally compute a
witness.

### Input

- `I1`      -- interval
- `I2`      -- interval
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``I1 ∩ I2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``I1 ∩ I2 = ∅``
  * `(false, v)` iff ``I1 ∩ I2 ≠ ∅`` and ``v ∈ I1 ∩ I2``

### Algorithm

``I1 ∩ I2 ≠ ∅`` iff there is a gap between the left-most point of the second
interval and the right-most point of the first interval, or vice-versa.

A witness is computed by taking the maximum over the left-most points of each
interval, which is guaranteed to belong to the intersection.
"""
function isdisjoint(I1::Interval, I2::Interval, witness::Bool=false)
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
    check = _isdisjoint(I1, I2, Val(false))
    if check
        N = promote_type(eltype(I1), eltype(I2))
        return (true, N[])
    else
        return (false, [max(min(I1), min(I2))])
    end
end
