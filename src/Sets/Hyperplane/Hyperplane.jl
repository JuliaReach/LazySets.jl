"""
    Hyperplane{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a hyperplane of the form ``a⋅x = b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The plane ``y = 0``:

```jldoctest
julia> Hyperplane([0, 1.], 0.)
Hyperplane{Float64, Vector{Float64}}([0.0, 1.0], 0.0)
```
"""
struct Hyperplane{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    function Hyperplane(a::VN, b::N) where {N,VN<:AbstractVector{N}}
        @assert !iszero(a) "a hyperplane needs a non-zero normal vector"
        return new{N,VN}(a, b)
    end
end
