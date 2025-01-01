export AbstractPolyhedron

"""
    AbstractPolyhedron{N} <: ConvexSet{N}

Abstract type for closed convex polyhedral sets.

### Notes

Every concrete `AbstractPolyhedron` must define the following functions:

- `constraints_list(::AbstractPolyhedron)` -- return a list of all facet constraints

Polyhedra are defined as the intersection of a finite number of closed
half-spaces.
As such, polyhedra are closed and convex but not necessarily bounded.
Bounded polyhedra are called *polytopes* (see [`AbstractPolytope`](@ref)).

The subtypes of `AbstractPolyhedron` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolyhedron)
8-element Vector{Any}:
 AbstractPolytope
 HPolyhedron
 HalfSpace
 Hyperplane
 Line
 Line2D
 Star
 Universe
```
"""
abstract type AbstractPolyhedron{N} <: ConvexSet{N} end

# supertype for all linear map algorithms
abstract type AbstractLinearMapAlgorithm end

"""
    constrained_dimensions(P::AbstractPolyhedron)

Return the indices in which a polyhedron is constrained.

### Input

- `P` -- polyhedron

### Output

A vector of ascending indices `i` such that the polyhedron is constrained in
dimension `i`.

### Examples

A 2D polyhedron with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(P::AbstractPolyhedron)
    constraints = constraints_list(P)
    if isempty(constraints)
        return Int[]
    end
    zero_indices = zeros(Int, dim(P))
    for constraint in constraints
        for i in constrained_dimensions(constraint)
            zero_indices[i] = i
        end
    end
    return filter(x -> x != 0, zero_indices)
end

# generic function for the AbstractPolyhedron interface => returns an HPolyhedron
function _linear_map_hrep_helper(M::AbstractMatrix, P::LazySet,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, P, algo)
    return HPolyhedron(constraints)
end

# internal functions; defined here due to dependency SymEngine and submodules
function isfeasible end
function _ishalfspace end
function _ishyperplanar end
function _parse_linear_expression end

# To account for the compilation order, other functions are defined in the file
# AbstractPolyhedron_functions.jl
