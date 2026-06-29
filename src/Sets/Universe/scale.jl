# see ext/LazySets/LazySetsUniverseExt.jl

function scale!(α::Real, U::Universe)
    if iszero(α)
        throw(ArgumentError("cannot 0-scale a universe in-place"))
    end
    return U
end
