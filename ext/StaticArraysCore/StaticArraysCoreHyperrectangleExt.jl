using LazySets: radius_hyperrectangle
using LazySets.HyperrectangleModule: Hyperrectangle
using ReachabilityBase.Comparison: isapproxzero
using StaticArraysCore: MMatrix, SMatrix, SVector
import LazySets: genmat, _genmat_static

function genmat(H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {L,N}
    gens = zeros(MMatrix{L,L,N})
    j = 0
    @inbounds for i in 1:L
        ri = radius_hyperrectangle(H, i)
        if !isapproxzero(ri)
            j += 1
            gens[i, j] = ri
        end
    end
    if j == L  # no zero radius
        return SMatrix{L,L}(gens)
    else
        return SMatrix{L,j}(view(gens, :, 1:j))::SMatrix{L}
    end
end

# this function is type stable, but it does not prune the generators
# according to flat dimensions of H
function _genmat_static(H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {L,N}
    gens = zeros(MMatrix{L,L,N})
    @inbounds for i in 1:L
        ri = radius_hyperrectangle(H, i)
        gens[i, i] = ri
    end
    return SMatrix{L,L}(gens)
end
