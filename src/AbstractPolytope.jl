import Base.isempty

export AbstractPolytope,
       vertices_list,
       singleton_list,
       linear_map,
       isempty

"""
    AbstractPolytope{N<:Real} <: LazySet{N}

Abstract type for polytopic sets, i.e., sets with finitely many flat facets, or
equivalently, sets defined as an intersection of a finite number of halfspaces,
or equivalently, sets with finitely many vertices.

### Notes

Every concrete `AbstractPolytope` must define the following functions:
- `constraints_list(::AbstractPolytope{N})::Vector{LinearConstraint{N}}` --
    return a list of all facet constraints
- `vertices_list(::AbstractPolytope{N})::Vector{Vector{N}}` -- return a list of
    all vertices

```jldoctest
julia> subtypes(AbstractPolytope)
4-element Array{Any,1}:
 AbstractCentrallySymmetricPolytope
 AbstractPolygon
 HPolytope
 VPolytope
```
"""
abstract type AbstractPolytope{N<:Real} <: LazySet{N} end


# --- common AbstractPolytope functions ---


"""
    singleton_list(P::AbstractPolytope{N})::Vector{Singleton{N}} where {N<:Real}

Return the vertices of a polytopic set as a list of singletons.

### Input

- `P` -- a polytopic set

### Output

List containing a singleton for each vertex.
"""
function singleton_list(P::AbstractPolytope{N}
                       )::Vector{Singleton{N}} where {N<:Real}
    return [Singleton(vi) for vi in vertices_list(P)]
end

"""
    linear_map(M::AbstractMatrix, P::AbstractPolytope{N};
               output_type::Type{<:LazySet}=VPolytope{N}) where {N<:Real}

Concrete linear map of an abstract polytype.

### Input

- `M`           -- matrix
- `P`           -- abstract polytype
- `output_type` -- (optional, default: `VPolytope`) type of the result

### Output

A set of type `output_type`.

### Algorithm

The linear map ``M`` is applied to each vertex of the given set ``P``, obtaining
a polytope in V-representation. Since some set representations (e.g. axis-aligned
hyperrectangles) are not closed under linear maps, the default output is a
`VPolytope`. If an `output_type` is given, the corresponding `convert` method
is invoked.
"""
function linear_map(M::AbstractMatrix, P::AbstractPolytope{N};
                    output_type::Type{<:LazySet}=VPolytope{N}) where {N<:Real}
    @assert dim(P) == size(M, 2)
    MP = broadcast(v -> M * v, vertices_list(P)) |> VPolytope{N}
    return convert(output_type, MP)
end

"""
    isempty(P::AbstractPolytope{N})::Bool where {N<:Real}

Determine whether a polytope is empty.

### Input

- `P` -- abstract polytope

### Output

`true` if the given polytope contains no vertices, and `false` otherwise.

### Algorithm

This algorithm checks whether the `vertices_list` of the given polytope is empty
or not.
"""
function isempty(P::AbstractPolytope{N})::Bool where {N<:Real}
    return isempty(vertices_list(P))
end

function default_polyhedra_backend(N)
    @assert isdefined(Main, :Polyhedra) "this function needs the package 'Polyhedra' to be loaded"
    error("no default backend for numeric type $N")
end

function load_polyhedra_abstractpolytope() # function to be loaded by Requires

return quote

import CDDLib # default backend
import Polyhedra:polyhedron,
                 HRep, VRep,
                 removehredundancy!, removevredundancy!,
                 hreps, vreps,
                 intersect,
                 convexhull,
                 hcartesianproduct,
                 points

@static if VERSION < v"0.7-"
    import Polyhedra:SimpleHRepresentation, SimpleHRepresentation
end

function default_polyhedra_backend(N::Type{<:AbstractFloat})
    return CDDLib.CDDLibrary()
end

function default_polyhedra_backend(N::Type{<:Rational})
    return CDDLib.CDDLibrary(:exact)
end

end # quote
end # function load_polyhedra_hpolytope()
