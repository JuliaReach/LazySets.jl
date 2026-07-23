import LazySets

using LazySets.Ball1Module: Ball1
using LazySets.Ball2Module: Ball2
using LazySets.BallInfModule: BallInf
import LazySets.BallpModule: _Ballp_special_cases

function _Ballp_special_cases(p::N, center::VN,
                              radius::N) where {N,VN<:AbstractVector{N}}
    if p == N(Inf)
        return BallInf(center, radius)
    elseif p == N(2)
        return Ball2(center, radius)
    elseif isone(p)
        return Ball1(center, radius)
    end
    return nothing
end
