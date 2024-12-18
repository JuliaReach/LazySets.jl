function ∈(x::AbstractVector, Z::ZeroSet)
    @assert length(x) == dim(Z) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(Z))-dimensional set"
    return iszero(x)
end
