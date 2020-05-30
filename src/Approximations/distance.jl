export distance

"""
    distance(H1::AbstractHyperrectangle{N}, H2::AbstractHyperrectangle{N};
             p::Real=2) where {N<:Real}

Compute the standard distance between two hyperrectangular sets.

### Input

- `H1` -- hyperrectangular set
- `H2` -- hyperrectangular set
- `p`  -- (optional; default: `2`) value of the ``p``-norm

### Output

The distance, which is zero if the sets intersect and otherwise the ``p``-norm
of the shortest line segment between any pair of points.
"""
function distance(H1::AbstractHyperrectangle{N},
                  H2::AbstractHyperrectangle{N};
                  p::Real=2) where {N<:Real}
    n = dim(H1)
    d = Vector{N}(undef, n)
    for i in 1:n
        # find relative position in dimension i
        # (if c1 == c2, the results are equivalent independent of the branch)
        if center(H1, i) >= center(H2, i)
            lhs = low(H1, i)
            rhs = high(H2, i)
        else
            lhs = low(H2, i)
            rhs = high(H1, i)
        end
        if lhs <= rhs
            d[i] = zero(N)
        else
            d[i] = rhs - lhs
        end
    end
    return norm(d, p)
end
