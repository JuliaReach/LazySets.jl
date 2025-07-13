@validate function ρ(d::AbstractVector, Z::ZeroSet)
    N = promote_type(eltype(d), eltype(Z))
    return zero(N)
end
