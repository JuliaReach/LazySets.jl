import LazySets

using LazySets: AbstractLinearMapAlgorithm, constrained_dimensions,
                isequivalent, _linear_map_hrep
using LazySets.HalfSpaceModule: _non_element_halfspace
using LazySets.HPolygonModule: HPolygon
using LazySets.HyperplaneModule: _constraints_list_hyperplane,
                                 _σ_hyperplane_halfspace
using LazySets.Line2DModule: Line2D
using LazySets.SingletonModule: Singleton
using LazySets.UniverseModule: Universe
import LazySets.API: constraints_list, isuniversal, project, σ
import LazySets: _linear_map_hrep_helper

function constraints_list(L::Line2D)
    return _constraints_list_hyperplane(L.a, L.b)
end

"""
# Extended help

    isuniversal(L::Line2D, [witness]::Bool=false)

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
function isuniversal(L::Line2D, witness::Bool=false)
    if witness
        v = _non_element_halfspace(L.a, L.b)
        return (false, v)
    else
        return false
    end
end

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Line2D{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a line  # TODO assert this
        c = first(constraints)
        return Line2D(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    elseif length(constraints) == 4
        # projection to the origin
        S = Singleton(zeros(N, 2))
        @assert isequivalent(HPolygon(constraints), S) "unexpected constraints"
        return S
    else
        throw(ArgumentError("unexpected number of $(length(constraints)) constraints"))
    end
end

# the algorithm is a 2D specialization of the `Hyperplane` algorithm, except
# that it returns a `Singleton` for a 1D line
@validate function project(L::Line2D{N}, block::AbstractVector{Int}; kwargs...) where {N}
    m = length(block)
    if m == 2
        @inbounds if block[1] == 1 && block[2] == 2
            return L  # no projection
        end
        # block[1] == 2 && block[2] == 1
        return Line2D(L.a[block], L.b)  # swap a vector
    end
    # projection to 1 dimension
    cdims = constrained_dimensions(L)
    if length(cdims) == 1
        @inbounds if cdims[1] == block[1]
            # L: aᵢxᵢ = b where aᵢ ≠ 0
            return Singleton([L.b / L.a[cdims[1]]])
        else
            # L: aⱼxⱼ = b where i ≠ j
            return Universe{N}(1)
        end
    end
    @assert length(cdims) == 2 "the line should be constrained in both dimensions"
    return Universe{N}(1)
end

@validate function σ(d::AbstractVector, L::Line2D)
    v, unbounded = _σ_hyperplane_halfspace(d, L.a, L.b; error_unbounded=true,
                                           halfspace=false)
    return v
end
