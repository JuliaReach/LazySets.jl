using LazySets: AbstractLinearMapAlgorithm, an_element, constrained_dimensions,
                _witness_result_empty, _linear_map_hrep, @validate
using LazySets.HalfSpaceModule: HalfSpace
using LazySets.HPolyhedronModule: HPolyhedron
using LazySets.HyperplaneModule: Hyperplane, _an_element_helper_hyperplane,
                                 _σ_hyperplane_halfspace
using LazySets.UniverseModule: Universe
using LinearAlgebra: dot
using ReachabilityBase.Arrays: samedir
import Base: isdisjoint
import LazySets: _linear_map_hrep_helper
import LazySets.API: project, ρ, σ
import LazySets.HalfSpaceModule: _an_element_halfspace

@inline function _an_element_halfspace(a::AbstractVector, b;
                                       nonzero_entry_a::Int=findfirst(!iszero, a),
                                       direction_inside::Bool)
    x = _an_element_helper_hyperplane(a, b, nonzero_entry_a)
    if direction_inside
        x .-= a
    else
        x .+= a
    end
    return x
end

function _linear_map_hrep_helper(M::AbstractMatrix, hs::HalfSpace,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, hs, algo)
    if length(constraints) == 1
        return first(constraints)
    elseif isempty(constraints)
        N = promote_type(eltype(M), eltype(hs))
        return Universe{N}(size(M, 1))
    else
        return HPolyhedron(constraints)
    end
end

"""
# Extended help

    project(H::HalfSpace, block::AbstractVector{Int}; [kwargs...])

### Algorithm

If the unconstrained dimensions of `H` are a subset of the `block` variables,
the projection is applied to the normal direction of `H`.
Otherwise, the projection results in the universal set.

The latter can be seen as follows.
Without loss of generality consider projecting out a single and constrained
dimension ``xₖ`` (projecting out multiple dimensions can be modeled by
repeatedly projecting out one dimension).
We can write the projection as an existentially quantified linear constraint:

```math
    ∃xₖ: a₁x₁ + … + aₖxₖ + … + aₙxₙ ≤ b
```

Since ``aₖ ≠ 0``, there is always a value for ``xₖ`` that satisfies the
constraint for any valuation of the other variables.

### Examples

Consider the half-space ``x + y + 0⋅z ≤ 1``, whose ambient dimension is `3`.
The (trivial) projection in the three dimensions using the block of variables
`[1, 2, 3]` is:

```jldoctest project_halfspace
julia> H = HalfSpace([1.0, 1.0, 0.0], 1.0)
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, 0.0], 1.0)

julia> project(H, [1, 2, 3])
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, 0.0], 1.0)
```

Projecting along dimensions `1` and `2` only:

```jldoctest project_halfspace
julia> project(H, [1, 2])
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```

For convenience, one can use `project(H, constrained_dimensions(H))` to return
the half-space projected on the dimensions where it is constrained:

```jldoctest project_halfspace
julia> project(H, constrained_dimensions(H))
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```

If a constrained dimension is projected, we get the universal set of the
dimension corresponding to the projection.

```jldoctest project_halfspace
julia> project(H, [1, 3])
Universe{Float64}(2)

julia> project(H, [1])
Universe{Float64}(1)
```
"""
@validate function project(H::HalfSpace, block::AbstractVector{Int}; kwargs...)
    if constrained_dimensions(H) ⊆ block
        return HalfSpace(H.a[block], H.b)
    else
        N = eltype(H)
        return Universe{N}(length(block))
    end
end

"""
# Extended help

    ρ(d::AbstractVector, hs::HalfSpace)

### Output

Unless the direction is (a multiple of) the normal direction of the half-space,
the result is `Inf`.
"""
@validate function ρ(d::AbstractVector, hs::HalfSpace)
    v, unbounded = _σ_hyperplane_halfspace(d, hs.a, hs.b; error_unbounded=false,
                                           halfspace=true)
    if unbounded
        N = promote_type(eltype(d), eltype(hs))
        return N(Inf)
    end
    return dot(d, v)
end

"""
# Extended help

    σ(d::AbstractVector, hs::HalfSpace)

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is (a multiple of) the normal direction of the half-space.
In both cases the result is any point on the boundary (the defining hyperplane).
Otherwise this function throws an error.
"""
@validate function σ(d::AbstractVector, hs::HalfSpace)
    v, unbounded = _σ_hyperplane_halfspace(d, hs.a, hs.b; error_unbounded=true,
                                           halfspace=true)
    return v
end

"""
# Extended help

    isdisjoint(H1::HalfSpace, H2::HalfSpace, [witness]::Bool=false)

### Algorithm

Two half-spaces do not intersect if and only if their normal vectors point in
the opposite direction and there is a gap between the two defining hyperplanes.

The latter can be checked as follows:
Let ``H1 : a_1⋅x = b_1`` and ``H2 : a_2⋅x = b_2``.
Then we already know that ``a_2 = -k⋅a_1`` for some positive scaling factor
``k``.
Let ``x_1`` be a point on the defining hyperplane of ``H1``.
We construct a line segment from ``x_1`` to the point ``x_2`` on the defining
hyperplane of ``hs_2`` by shooting a ray from ``x_1`` with direction ``a_1``.
Thus we look for a factor ``s`` such that ``(x_1 + s⋅a_1)⋅a_2 = b_2``.
This gives us ``s = (b_2 - x_1⋅a_2) / (-k a_1⋅a_1)``.
The gap exists if and only if ``s`` is positive.

If the normal vectors do not point in opposite directions, then the defining
hyperplanes intersect and we can produce a witness as follows.
All points ``x`` in this intersection satisfy ``a_1⋅x = b_1`` and
``a_2⋅x = b_2``. Thus we have ``(a_1 + a_2)⋅x = b_1+b_2``.
We now find a dimension where ``a_1 + a_2`` is non-zero, say, ``i``.
Then the result is a vector with one non-zero entry in dimension ``i``, defined
as ``[0, …, 0, (b_1 + b_2)/(a_1[i] + a_2[i]), 0, …, 0]``.
Such a dimension ``i`` always exists.
"""
@validate function isdisjoint(H1::HalfSpace, H2::HalfSpace, witness::Bool=false)
    a1 = H1.a
    a2 = H2.a
    N = promote_type(eltype(H1), eltype(H2))
    issamedir, k = samedir(a1, -a2)
    if issamedir
        x1 = an_element(Hyperplane(a1, H1.b))
        b2 = H2.b
        s = (b2 - dot(x1, a2)) / (-k * dot(a1, a1))
        empty_intersection = s > 0
        # if `!empty_intersection`, x1 is a witness because both defining
        # hyperplanes are contained in each half-space
        return _witness_result_empty(witness, empty_intersection, H1, H2, x1)
    elseif !witness
        return false
    end
    # compute witness
    v = zeros(N, length(a1))
    for i in eachindex(a1)
        a_sum_i = a1[i] + a2[i]
        if !iszero(a_sum_i)
            v[i] = (H1.b + H2.b) / a_sum_i
            break
        end
    end
    return (false, v)
end
