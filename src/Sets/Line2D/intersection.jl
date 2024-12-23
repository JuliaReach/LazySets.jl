"""
# Extended help

    intersection(L1::Line2D, L2::Line2D)

### Output

Three outcomes are possible:

- If the lines are identical, the result is the first line.
- If the lines are parallel and not identical, the result is the empty set.
- Otherwise the result is the set with the unique intersection point.

### Algorithm

We first check whether the lines are parallel.
If not, we use [Cramer's rule](https://en.wikipedia.org/wiki/Cramer%27s_rule)
to compute the intersection point.

### Examples

The line ``y = x`` intersected with the line ``y = -x + 1`` respectively with
itself:

```jldoctest
julia> intersection(Line2D([-1.0, 1], 0.0), Line2D([1.0, 1], 1.0))
Singleton{Float64, Vector{Float64}}([0.5, 0.5])

julia> intersection(Line2D([1.0, 1], 1.0), Line2D([1.0, 1], 1.0))
Line2D{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```
"""
function intersection(L1::Line2D, L2::Line2D)
    return _intersection_line2d(L1, L2)
end
