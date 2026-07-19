import LazySets

using LazySets: @validate
using LazySets.HalfSpaceModule: HalfSpace
using LazySets.LineModule: Line
using LazySets.SingletonModule: Singleton
using LazySets.UniverseModule: Universe
using LinearAlgebra: nullspace
import LazySets.API: constraints_list, project, linear_map

"""
# Extended help

    constraints_list(L::Line)

### Output

A list containing `2n-2` half-spaces whose intersection is `L`, where `n` is the
ambient dimension of `L`.
"""
function constraints_list(L::Line)
    p = L.p
    n = length(p)
    d = reshape(L.d, 1, n)
    K = nullspace(d)
    m = size(K, 2)
    @assert m == n - 1 "expected $(n - 1) normal half-spaces, got $m"

    N, VN = _parameters(L)
    out = Vector{HalfSpace{N,VN}}(undef, 2m)
    idx = 1
    @inbounds for Kj in eachcol(K)
        b = dot(Kj, p)
        out[idx] = HalfSpace(Kj, b)
        out[idx + 1] = HalfSpace(-Kj, -b)
        idx += 2
    end
    return out
end

# reason: Documenter cannot annotate `constraints_list` with type parameters
function _parameters(::Line{N,VN}) where {N,VN}
    return (N, VN)
end

"""
# Extended help

    linear_map(M::AbstractMatrix, L::Line)

### Output

The line obtained by applying the linear map, if that still results in a line,
or a `Singleton` otherwise.

### Algorithm

We apply the linear map to the point and direction of `L`.
If the resulting direction is zero, the result is a singleton.
"""
@validate function linear_map(M::AbstractMatrix, L::Line)
    Mp = M * L.p
    Md = M * L.d
    if iszero(Md)
        return Singleton(Mp)
    end
    return Line(Mp, Md)
end

@validate function project(L::Line{N}, block::AbstractVector{Int}; kwargs...) where {N}
    d = L.d[block]
    if iszero(d)
        return Singleton(L.p[block])  # projected out all nontrivial dimensions
    elseif length(d) == 1
        return Universe{N}(1)  # special case: 1D line is a universe
    else
        return Line(L.p[block], d)
    end
end
