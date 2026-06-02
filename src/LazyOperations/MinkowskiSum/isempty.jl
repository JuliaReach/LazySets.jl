"""
    isempty(ms::MinkowskiSum)

Check whether a Minkowski sum of two sets is empty.

### Input

- `ms` -- Minkowski sum of two sets

### Output

`true` iff any of the wrapped sets are empty.
"""
function isempty(ms::MinkowskiSum)
    return isempty(ms.X) || isempty(ms.Y)
end
