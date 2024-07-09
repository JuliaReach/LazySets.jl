"""
    ZeroSet{N} <: AbstractSingleton{N}

Type that represents the zero set, i.e., the set that only contains the origin.

### Fields

- `dim` -- the ambient dimension of the set
"""
struct ZeroSet{N} <: AbstractSingleton{N}
    dim::Int
end

# default constructor of type Float64
ZeroSet(dim::Int) = ZeroSet{Float64}(dim)
