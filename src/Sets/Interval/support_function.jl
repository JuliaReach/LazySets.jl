function ρ(d::AbstractVector, x::Interval)
    @assert length(d) == dim(x) "a $(length(d))-dimensional vector is " *
                                "incompatible with an $(dim(x))-dimensional set"
    N = promote_type(eltype(d), eltype(x))
    return @inbounds d[1] * (d[1] > zero(N) ? max(x) : min(x))
end
