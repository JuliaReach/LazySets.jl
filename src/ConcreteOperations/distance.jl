export distance

"""
    distance(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle;
             [p]::Real=2)

Compute the standard distance between two hyperrectangular sets, defined as

```math
    \\inf_{x \\in H_1, y \\in H_2} \\{ d(x, y) \\}.
```

### Input

- `H1` -- hyperrectangular set
- `H2` -- hyperrectangular set
- `p`  -- (optional; default: `2`) value of the ``p``-norm

### Output

The distance, which is zero if the sets intersect and otherwise the ``p``-norm
of the shortest line segment between any pair of points.

### Notes

See also [`hausdorff_distance`](@ref) for an alternative distance notion.
"""
function distance(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle;
                  p::Real=2)
    n = dim(H1)
    @assert n == dim(H2) "incompatible set dimensions $n and $(dim(H2))"

    N = promote_type(eltype(H1), eltype(H2))
    d = Vector{N}(undef, n)
    @inbounds for i in 1:n
        # find relative position in dimension i
        # (if c1 == c2, the results are equivalent independent of the branch)
        if center(H1, i) >= center(H2, i)
            lhs = low(H1, i)
            rhs = high(H2, i)
        else
            lhs = low(H2, i)
            rhs = high(H1, i)
        end
        if _leq(lhs, rhs)
            d[i] = zero(N)
        else
            d[i] = rhs - lhs
        end
    end
    return norm(d, p)
end

"""
    distance(S::AbstractSingleton, X::ConvexSet; [p]::Real=2.0)

Compute the distance between the singleton `S` and the set `X` with respect to
the given `p`-norm.

### Input

- `S` -- singleton, i.e., a set with one element
- `X` -- set
- `p` -- (optional, default: `2.0`) the `p`-norm used; `p = 2.0` corresponds to
         the usual Euclidean norm

### Output

A scalar representing the distance between the element wrapped by `S` and the
set `X`.
"""
function distance(S::AbstractSingleton, X::ConvexSet; p::Real=2.0)
    return distance(element(S), X; p=p)
end

distance(X::ConvexSet, S::AbstractSingleton; p::Real=2.0) = distance(S, X; p=p)

distance(S1::AbstractSingleton, S2::AbstractSingleton; p::Real=2.0) =
    distance(element(S1), element(S2); p=p)

# disambiguation
distance(S::AbstractSingleton, H::AbstractHyperrectangle; p::Real=2.0) =
    distance(element(S), H; p=p)
distance(H::AbstractHyperrectangle, S::AbstractSingleton; p::Real=2.0) =
    distance(element(S), H; p=p)
