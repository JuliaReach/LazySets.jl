function ∈(x::AbstractVector, U::Universe)
    @assert length(x) == dim(U) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(U))-dimensional set"
    return true
end
