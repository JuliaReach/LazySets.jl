"""
# Extended help

    ∈(x::AbstractVector, B::Ball1, [failfast]::Bool=false)

### Input

- `x` -- point/vector
- `B` -- ball in the 1-norm
- `failfast` -- (optional, default: `false`) optimization for negative answer

### Notes

The default behavior (`failfast == false`) is worst-case optimized, i.e., the
implementation is optimistic and first computes (see below) the whole sum before
comparing to the radius.
In applications where the point is typically far away from the ball, the option
`failfast == true` terminates faster.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the 1-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``∑_{i=1}^n |c_i - x_i| ≤ r``.

### Examples

```jldoctest
julia> B = Ball1([1.0, 1.0], 1.0);

julia> [0.5, -0.5] ∈ B
false
julia> [0.5, 1.5] ∈ B
true
```
"""
function ∈(x::AbstractVector, B::Ball1, failfast::Bool=false)
    @assert length(x) == dim(B) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    N = promote_type(eltype(x), eltype(B))
    sum = zero(N)
    @inbounds for (i, xi) in enumerate(x)
        sum += abs(B.center[i] - xi)
        if failfast && sum > B.radius
            return false
        end
    end
    return sum <= B.radius
end
