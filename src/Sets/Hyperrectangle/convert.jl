function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

function convert(::Type{Hyperrectangle}, I::IA.Interval)
    low_I = [IA.inf(I)]
    high_I = [IA.sup(I)]
    return Hyperrectangle(; low=low_I, high=high_I)
end
