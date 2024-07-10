function scale(α::Real, U::Universe{N}) where {N}
    if iszero(α)
        require(@__MODULE__, :LazySets; fun_name="scale")

        return ZeroSet{N}(dim(U))
    end
    return U
end

function scale!(α::Real, U::Universe)
    if iszero(α)
        throw(ArgumentError("cannot 0-scale a universe in-place"))
    end
    return U
end
