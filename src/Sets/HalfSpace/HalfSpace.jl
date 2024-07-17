"""
    HalfSpace{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a (closed) half-space of the form ``a⋅x ≤ b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The half-space ``x + 2y - z ≤ 3``:

```jldoctest
julia> HalfSpace([1, 2, -1.], 3.)
HalfSpace{Float64, Vector{Float64}}([1.0, 2.0, -1.0], 3.0)
```

To represent the set ``y ≥ 0`` in the plane, we can rearrange the expression as
``0x - y ≤ 0``:

```jldoctest
julia> HalfSpace([0, -1.], 0.)
HalfSpace{Float64, Vector{Float64}}([0.0, -1.0], 0.0)
```
"""
struct HalfSpace{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    function HalfSpace(a::VN, b::N) where {N,VN<:AbstractVector{N}}
        @assert !iszero(a) "a half-space needs a non-zero normal vector"
        return new{N,VN}(a, b)
    end
end
