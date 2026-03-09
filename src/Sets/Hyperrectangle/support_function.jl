@validate function ρ(d::SingleEntryVector, H::Hyperrectangle)
    return _ρ_sev_hyperrectangle(d, H)
end

# hyperrectangle radius given as single-entry vector
# (i.e., hyperrectangle is flat in all but one direction)
@validate function ρ(d::AbstractVector, H::Hyperrectangle{N,VNC,<:SingleEntryVector}) where {N,VNC}
    N2 = promote_type(eltype(d), eltype(H))
    r = radius_hyperrectangle(H)
    res = dot(d, center(H))
    di = @inbounds d[r.i]
    if di < zero(N2)
        res -= di * r.v
    elseif di > zero(N2)
        res += di * r.v
    end
    return res
end

# both direction and hyperrectangle radius given as single-entry vector
@validate function ρ(d::SingleEntryVector,
                     H::Hyperrectangle{N,VNC,SingleEntryVector{N}}) where {N,VNC}
    N2 = promote_type(eltype(d), eltype(H))
    r = radius_hyperrectangle(H)
    res = d.v * center(H, d.i)
    if d.i == r.i
        if d.v < zero(N2)
            res -= d.v * r.v
        elseif d.v > zero(N2)
            res += d.v * r.v
        end
    end
    return res
end
