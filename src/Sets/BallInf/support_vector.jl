@validate function σ(d::AbstractVector, B::BallInf)
    return center(B) .+ sign.(d) .* B.radius
end

# special case for SingleEntryVector
@validate function σ(d::SingleEntryVector, B::BallInf)
    return _σ_sev_hyperrectangle(d, B)
end
