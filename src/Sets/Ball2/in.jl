"""
# Extended help

    ∈(x::AbstractVector, B::Ball2)

### Notes

This implementation is worst-case optimized, i.e., it is optimistic and first
computes (see below) the whole sum before comparing to the radius.
In applications where the point is typically far away from the ball, a fail-fast
implementation with interleaved comparisons could be more efficient.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the 2-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``\\left( ∑_{i=1}^n |c_i - x_i|^2 \\right)^{1/2} ≤ r``.

### Examples

```jldoctest
julia> B = Ball2([1., 1.], sqrt(0.5))
Ball2{Float64, Vector{Float64}}([1.0, 1.0], 0.7071067811865476)

julia> [.5, 1.6] ∈ B
false

julia> [.5, 1.5] ∈ B
true
```
"""
function ∈(x::AbstractVector, B::Ball2)
    @assert length(x) == dim(B)
    N = promote_type(eltype(x), eltype(B))
    sum = zero(N)
    @inbounds for i in eachindex(x)
        sum += (B.center[i] - x[i])^2
    end
    return _leq(sqrt(sum), B.radius)
end
