import LazySets

using LazySets: AbstractLinearMapAlgorithm, an_element, constrained_dimensions,
                intersection, isequivalent, _linear_map_hrep,
                _witness_result_empty, @validate_commutative
using LazySets.HalfSpaceModule: HalfSpace, _non_element_halfspace,
                                _normalize_halfspace
using LazySets.HPolyhedronModule: HPolyhedron
using LazySets.HyperplaneModule: Hyperplane
using LazySets.EmptySetModule: EmptySet
using LazySets.UniverseModule: Universe
using LinearAlgebra: dot
import LazySets: normalize, _linear_map_hrep_helper
import LazySets.API: distance, isuniversal, project
import LazySets.HyperplaneModule: _constraints_list_hyperplane,
                                  _isdisjoint_hyperplane_hyperplane

function _constraints_list_hyperplane(a::AbstractVector, b)
    return [HalfSpace(a, b), HalfSpace(-a, -b)]
end

"""
# Extended help

    isuniversal(H::Hyperplane, [witness]::Bool=false)

### Algorithm

A witness is produced by adding the normal vector to an element on the
hyperplane.
"""
function isuniversal(H::Hyperplane, witness::Bool=false)
    if witness
        v = _non_element_halfspace(H.a, H.b)
        return (false, v)
    else
        return false
    end
end

@validate_commutative function distance(x::AbstractVector, H::Hyperplane; p::Real=2)
    if p != 2
        throw(ArgumentError("`distance` is only implemented for Euclidean norm"))
    end

    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return abs(dot(x, a) - b)
end

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Hyperplane{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a hyperplane
        c = first(constraints)
        return Hyperplane(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    else
        return HPolyhedron(constraints)
    end
end

@validate function project(H::Hyperplane{N}, block::AbstractVector{Int}; kwargs...) where {N}
    if constrained_dimensions(H) ⊆ block
        return Hyperplane(H.a[block], H.b)
    else
        return Universe{N}(length(block))
    end
end

function _isdisjoint_hyperplane_hyperplane(hp1::Hyperplane,
                                           hp2::Hyperplane,
                                           witness::Bool=false)
    if isequivalent(hp1, hp2)
        res = false
        if witness
            w = an_element(hp1)
        end
    else
        cap = intersection(hp1, hp2)
        res = cap isa EmptySet
        if !res && witness
            w = an_element(cap)
        end
    end
    if res
        return _witness_result_empty(witness, true, hp1, hp2)
    end
    return witness ? (false, w) : false
end

"""
    normalize(H::Hyperplane{N}, p::Real=N(2)) where {N}

Normalize a hyperplane.

### Input

- `H` -- hyperplane
- `p` -- (optional, default: `2`) norm

### Output

A new hyperplane whose normal direction ``a`` is normalized, i.e., such that
``‖a‖_p = 1`` holds.
"""
function normalize(H::Hyperplane{N}, p::Real=N(2)) where {N}
    a, b = _normalize_halfspace(H, p)
    return Hyperplane(a, b)
end
