@validate function in(x::AbstractVector, Z::ZeroSet)
    return iszero(x)
end
