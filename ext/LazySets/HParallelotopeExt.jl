import LazySets

using LazySets.HalfSpaceModule: HalfSpace
using LazySets.HPolyhedronModule: HPolyhedron
import LazySets.HParallelotopeModule: _constraints_list_hparallelotope,
                                      _HPolyhedron

function _constraints_list_hparallelotope(D, c::AbstractVector, N, VN)
    if isempty(D)
        return Vector{HalfSpace{N,VN}}(undef, 0)
    end
    n = size(D, 1)
    clist = Vector{HalfSpace{N,VN}}(undef, 2n)
    @inbounds for i in 1:n
        clist[i] = HalfSpace(D[i, :], c[i])
        clist[i + n] = HalfSpace(-D[i, :], c[i + n])
    end
    return clist
end

function _HPolyhedron(clist::AbstractVector)
    return HPolyhedron(clist)
end
