"""
    isbounded(ms::MinkowskiSum)

Check whether a Minkowski sum of two sets is bounded.

### Input

- `ms` -- Minkowski sum of two sets

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(ms::MinkowskiSum)
    return isbounded(ms.X) && isbounded(ms.Y)
end
