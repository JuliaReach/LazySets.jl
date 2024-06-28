function ∈(v::AbstractVector, x::Interval)
    @assert length(v) == 1 "a $(length(v))-dimensional vector is incompatible with an interval"
    return @inbounds v[1] ∈ x.dat
end
