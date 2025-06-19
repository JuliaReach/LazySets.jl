"""
# Extended help

    an_element(X::Interval)

### Algorithm

Return the left border (`low(X)`) of the interval.
"""
function an_element(X::Interval)
    return low(X)
end
