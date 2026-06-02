"""
    isbounded(em::ExponentialMap)

Check whether an exponential map is bounded.

### Input

- `em` -- exponential map

### Output

`true` iff the exponential map is bounded.
"""
function isbounded(em::ExponentialMap)
    return isbounded(em.X)
end
