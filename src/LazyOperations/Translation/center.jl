"""
    center(tr::Translation)

Return the center of the translation of a centrally-symmetric set.

### Input

- `tr` -- translation of a centrally-symmetric set

### Output

The translation of the center of the wrapped set by the translation vector.
"""
function center(tr::Translation)
    return center(tr.X) + tr.v
end
