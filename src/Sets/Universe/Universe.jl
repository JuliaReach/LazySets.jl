"""
    Universe{N} <: AbstractPolyhedron{N}

Type that represents the universal set, i.e., the set of all elements.

### Fields

- `dim` -- the ambient dimension of the set
"""
struct Universe{N} <: AbstractPolyhedron{N}
    dim::Int
end

# default constructor of type Float64
Universe(dim::Int) = Universe{Float64}(dim)
