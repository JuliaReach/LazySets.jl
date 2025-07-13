@validate function ∈(x::AbstractVector, Z::ZeroSet)
    return iszero(x)
end
