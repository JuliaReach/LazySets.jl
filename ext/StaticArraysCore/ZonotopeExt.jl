using StaticArraysCore: MMatrix, MVector, SMatrix, SVector
using LazySets: AbstractReductionMethod, GIR05, genmat
using LazySets.ZonotopeModule: Zonotope
import LazySets: reduce_order
import LazySets.ZonotopeModule: _split_ret

# conversion for static matrix
function reduce_order(Z::Zonotope{N,<:AbstractVector,<:MMatrix}, r::Real,
                      method::AbstractReductionMethod=GIR05()) where {N}
    return reduce_order(Zonotope(center(Z), SMatrix(genmat(Z))), r, method)
end

function _split_ret(Z₁::Zonotope{N,SV,SM},
                    Z₂::Zonotope{N,SV,SM}) where {N,n,p,SV<:MVector{n,N},SM<:MMatrix{n,p,N}}
    Z₁ = Zonotope(SVector{n}(Z₁.center), SMatrix{n,p}(Z₁.generators))
    Z₂ = Zonotope(SVector{n}(Z₂.center), SMatrix{n,p}(Z₂.generators))
    return Z₁, Z₂
end
