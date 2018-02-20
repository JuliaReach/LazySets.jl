#= concrete implementations of binary intersections between sets =#

export intersection

"""
    intersection(L1::Line{N}, L2::Line{N})::Vector{N} where {N<:Real}

Return the intersection of two 2D lines.

### Input

- `L1` -- first line

- `L2` -- second line

### Output

If the lines are parallel or identical, the result is an empty vector.
Otherwise the result is the only intersection point.

### Examples

The line ``y = -x + 1`` intersected with the line ``y = x``:

```jldoctest
julia> intersection(Line([-1., 1.], 0.), Line([1., 1.], 1.))
2-element Array{Float64,1}:
 0.5
 0.5
julia> intersection(Line([1., 1.], 1.), Line([1., 1.], 1.))
0-element Array{Float64,1}

```
"""
function intersection(L1::Line{N}, L2::Line{N})::Vector{N} where {N<:Real}
    b = [L1.b, L2.b]
    a = [L1.a.'; L2.a.']
    try
        # results in LAPACKException or SingularException if parallel
        return a \ b
    catch
        return N[]
    end
end
