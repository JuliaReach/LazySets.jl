"""
    Singleton{N, VN<:AbstractVector{N}} <: AbstractSingleton{N}

Type that represents a singleton, that is, a set with a single element.

### Fields

- `element` -- the only element of the set

### Examples

```jldoctest
julia> Singleton([1.0, 2.0])
Singleton{Float64, Vector{Float64}}([1.0, 2.0])

julia> Singleton(1.0, 2.0)  # convenience constructor from numbers
Singleton{Float64, Vector{Float64}}([1.0, 2.0])
```
"""
struct Singleton{N,VN<:AbstractVector{N}} <: AbstractSingleton{N}
    element::VN
end

function Singleton(args::Real...)
    return Singleton(collect(args))
end
