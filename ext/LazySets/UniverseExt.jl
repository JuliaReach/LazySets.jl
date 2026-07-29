using LazySets: LazySet, dim, isuniversal
using LazySets.EmptySetModule: EmptySet
using LazySets.HalfSpaceModule: HalfSpace
using LazySets.HPolyhedronModule: HPolyhedron
using LazySets.UniverseModule: Universe
using LazySets.ZeroSetModule: ZeroSet
using ReachabilityBase.Arrays: SingleEntryVector
import LazySets.API: complement, constraints_list, rectify, scale
import LazySets.UniverseModule: _difference_universe2,
                                _minkowski_difference_universe2

function complement(U::Universe{N}) where {N}
    return EmptySet{N}(dim(U))
end

function constraints_list(::Universe{N}) where {N}
    return HalfSpace{N,Vector{N}}[]
end

function rectify(U::Universe)
    N = eltype(U)
    n = dim(U)
    clist = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, n)
    @inbounds for i in 1:n
        clist[i] = HalfSpace(SingleEntryVector(i, n, N(-1)), zero(N))
    end
    return HPolyhedron(clist)
end

function scale(α::Real, U::Universe{N}) where {N}
    if iszero(α)
        return ZeroSet{N}(dim(U))
    end
    return U
end

function _difference_universe2(X::LazySet, U::Universe)
    N = promote_type(eltype(X), eltype(U))
    return EmptySet{N}(dim(U))
end

function _minkowski_difference_universe2(X::LazySet, U::Universe)
    if isuniversal(X)
        return U
    end
    N = promote_type(eltype(X), eltype(U))
    return EmptySet{N}(dim(U))
end
