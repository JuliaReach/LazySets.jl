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
    a = [transpose(L1.a); transpose(L2.a)]
    try
        # results in LAPACKException or SingularException if parallel
        return a \ b
    catch
        return N[]
    end
end

"""
    intersection(H1::AbstractHyperrectangle{N},
                 H2::AbstractHyperrectangle{N}
                )::Union{<:Hyperrectangle{N}, EmptySet{N}} where {N<:Real}

Return the intersection of two hyperrectangles.

### Input

- `H1` -- first hyperrectangle
- `H2` -- second hyperrectangle

### Output

If the hyperrectangles do not intersect, the result is the empty set.
Otherwise the result is the hyperrectangle that describes the intersection.

### Algorithm

In each isolated direction `i` we compute the rightmost left border and the
leftmost right border of the hyperrectangles.
If these borders contradict, then the intersection is empty.
Otherwise the result uses these borders in each dimension.
"""
function intersection(H1::AbstractHyperrectangle{N},
                      H2::AbstractHyperrectangle{N}
                     )::Union{<:Hyperrectangle{N}, EmptySet{N}} where {N<:Real}
    n = dim(H1)
    c1 = center(H1)
    c2 = center(H2)
    r1 = radius_hyperrectangle(H1)
    r2 = radius_hyperrectangle(H2)
    high = Vector{N}(n)
    low = Vector{N}(n)
    for i in 1:n
        high1 = c1[i] + r1[i]
        low1 = c1[i] - r1[i]
        high2 = c2[i] + r2[i]
        low2 = c2[i] - r2[i]
        high[i] = min(high1, high2)
        low[i] = max(low1, low2)
        if high[i] < low[i]
            return EmptySet{N}()
        end
    end
    return Hyperrectangle(high=high, low=low)
end
