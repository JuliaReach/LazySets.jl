"""
# Extended help

    an_element(x::Interval)

### Algorithm

Return the left border (`low(x)`) of the interval.
"""
function an_element(x::Interval)
    return low(x)
end
