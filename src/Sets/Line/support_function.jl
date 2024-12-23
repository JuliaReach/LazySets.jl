function ρ(d::AbstractVector, L::Line)
    if isapproxzero(dot(d, L.d))
        return dot(d, L.p)
    else
        N = eltype(L)
        return N(Inf)
    end
end
