function ρ(d::AbstractVector, X::Interval)
    @assert length(d) == dim(X) "a $(length(d))-dimensional vector is " *
                                "incompatible with an $(dim(X))-dimensional set"
    N = promote_type(eltype(d), eltype(X))
    return @inbounds d[1] * (d[1] > zero(N) ? max(X) : min(X))
end
