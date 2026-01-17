export AbstractPolytope,
       remove_redundant_vertices,
       remove_redundant_vertices!

"""
    AbstractPolytope{N} <: AbstractPolyhedron{N}

Abstract type for compact convex polytopic sets.

### Notes

See [`HPolytope`](@ref) or [`VPolytope`](@ref) for standard implementations of
this interface.

Every concrete `AbstractPolytope` must define the following method:

- `vertices_list(::AbstractPolytope)` -- return a list of all vertices

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolytope)
5-element Vector{Any}:
 AbstractCentrallySymmetricPolytope
 AbstractPolygon
 HPolytope
 Tetrahedron
 VPolytope
```

A polytope is a bounded polyhedron (see [`AbstractPolyhedron`](@ref)).
Polytopes are compact convex sets with either of the following equivalent
properties:
1. They are the intersection of a finite number of closed half-spaces.
2. They are the convex hull of finitely many vertices.
"""
abstract type AbstractPolytope{N} <: AbstractPolyhedron{N} end

function isboundedtype(::Type{<:AbstractPolytope})
    return true
end

function isbounded(::AbstractPolytope)
    return true
end

"""
# Extended help

    isempty(P::AbstractPolytope)

### Algorithm

This algorithm checks whether the `vertices_list` of `P` is empty.
"""
function isempty(P::AbstractPolytope)
    return isempty(vertices_list(P))
end

"""
# Extended help

    isuniversal(P::AbstractPolytope, [witness]::Bool=false)

### Algorithm

A witness is produced using `isuniversal(H)` where `H` is the first linear
constraint of `P`.
"""
function isuniversal(P::AbstractPolytope, witness::Bool=false)
    if witness
        constraints = constraints_list(P)
        if isempty(constraints)
            # special case for polytopes without constraints
            throw(ArgumentError("illegal polytope without constraints"))
        end
        return isuniversal(constraints[1], true)
    else
        return false
    end
end

# given a polytope P, apply the linear map P to each vertex of P
# it is assumed that the interface function `vertices_list(P)` is available
@inline function _linear_map_vrep(M::AbstractMatrix, P::AbstractPolytope,
                                  algo::LinearMapVRep=LinearMapVRep(nothing);
                                  apply_convex_hull::Bool=false)
    vlist = broadcast(v -> M * v, vertices_list(P))

    m = size(M, 1) # output dimension
    if m == 1
        convex_hull!(vlist)
        # points are sorted
        return Interval(vlist[1][1], vlist[end][1])
    elseif m == 2
        return VPolygon(vlist)
    else
        if apply_convex_hull
            convex_hull!(vlist)
        end
        return VPolytope(vlist)
    end
end

function _linear_map_hrep_helper(M::AbstractMatrix, P::AbstractPolytope,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, P, algo)
    m = size(M, 1) # output dimension
    if m == 1
        # TODO: create interval directly ?
        return convert(Interval, HPolytope(constraints))
    elseif m == 2
        return HPolygon(constraints)
    else
        return HPolytope(constraints)
    end
end

"""
# Extended help

    volume(P::AbstractPolytope; backend=default_polyhedra_backend(P))

### Input

- `backend` -- (optional, default: `default_polyhedra_backend(P)`) the backend
               for polyhedral computations; see [Polyhedra's
               documentation](https://juliapolyhedra.github.io/) for further
               information

### Algorithm

The volume is computed by the `Polyhedra` library.
"""
function volume(P::AbstractPolytope; backend=nothing)
    n = dim(P)
    if n <= 2
        if n == 1
            return _volume_1D(P)
        elseif n == 2
            return area(P)
        end
        throw(ArgumentError("invalid dimension $n"))
    end

    require(@__MODULE__, :Polyhedra; fun_name="volume")
    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    return Polyhedra.volume(polyhedron(P; backend=backend))
end

"""
    remove_redundant_vertices(P::AbstractPolytope)

Return an equivalent polytope in vertex representation with redundant vertices
removed.

### Input

- `P` -- polytope in vertex representation

### Output

A new polytope with the redundant vertices removed.
"""
function remove_redundant_vertices(::AbstractPolytope) end

"""
    remove_redundant_vertices!(P::AbstractPolytope)

Remove the redundant vertices from a polytope in vertex representation in-place.

### Input

- `P` -- polytope in vertex representation

### Output

A new polytope with the redundant vertices removed.
"""
function remove_redundant_vertices!(::AbstractPolytope) end

function constrained_dimensions(P::AbstractPolytope)
    return 1:dim(P)
end
