module StaticArraysCoreExt

using LazySets: center, _genmat_static
using LazySets.HyperrectangleModule: Hyperrectangle
using LazySets.ZonotopeModule: Zonotope
using StaticArraysCore: MMatrix, SMatrix, SVector
import LazySets: _convert_static, _hcat_KLred, _interval_hull, _to_colVector
import LazySets.Approximations: dir_east, dir_north, dir_west, dir_south

include("StaticArraysCore/StaticArraysCoreBallInfExt.jl")
include("StaticArraysCore/StaticArraysCoreHyperrectangleExt.jl")
include("StaticArraysCore/StaticArraysCoreZonotopeExt.jl")

const DIR_EAST_STATIC(N) = SVector{2}([one(N), zero(N)])
const DIR_NORTH_STATIC(N) = SVector{2}([zero(N), one(N)])
const DIR_WEST_STATIC(N) = SVector{2}([-one(N), zero(N)])
const DIR_SOUTH_STATIC(N) = SVector{2}([zero(N), -one(N)])

dir_east(N, ::SVector) = DIR_EAST_STATIC(N)
dir_north(N, ::SVector) = DIR_NORTH_STATIC(N)
dir_west(N, ::SVector) = DIR_WEST_STATIC(N)
dir_south(N, ::SVector) = DIR_SOUTH_STATIC(N)

function _to_colVector(M::SMatrix)
    return collect(eachcol(M))
end

function _convert_static(::Type{Zonotope}, H::Hyperrectangle{N,<:SVector,<:SVector}) where {N}
    return Zonotope(center(H), _genmat_static(H))
end

# implementation for static arrays
function _interval_hull(G::SMatrix{n,p,N,L}, indices) where {n,p,N,L}
    Lred = zeros(MMatrix{n,n,N})
    @inbounds for i in 1:n
        for j in indices
            Lred[i, i] += abs(G[i, j])
        end
    end
    return SMatrix{n,n}(Lred)
end

# implementation for static arrays
function _hcat_KLred(G::SMatrix{n,p,N,L1}, indices, Lred::SMatrix{n,n,N,L2}) where {n,p,N,L1,L2}
    m = length(indices)
    K = SMatrix{n,m}(view(G, :, indices))
    return hcat(K, Lred)
end

end  # module
