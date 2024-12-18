function σ(d::AbstractVector, B::BallInf)
    @assert length(d) == dim(B) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    return center(B) .+ sign.(d) .* B.radius
end

# special case for SingleEntryVector
function σ(d::SingleEntryVector, B::BallInf)
    return _σ_sev_hyperrectangle(d, B)
end
