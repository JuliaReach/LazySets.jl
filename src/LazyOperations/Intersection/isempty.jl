"""
    isempty(cap::Intersection)

Check whether the intersection of two sets is empty.

### Input

- `cap` -- intersection of two sets

### Output

`true` iff the intersection is empty.

### Notes

The result will be cached, so a second query will be fast.
"""
function isempty(cap::Intersection)
    if isempty_known(cap)
        # use cached result
        return isempty(cap.cache)
    end
    # compute result
    empty_intersection = isdisjoint(cap.X, cap.Y)
    # update cache
    set_isempty!(cap, empty_intersection)

    return empty_intersection
end
