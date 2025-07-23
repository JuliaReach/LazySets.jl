function scale(α::Real, MZ::MatrixZonotope)
    return scale!(α, copy(MZ))
end

function scale!(α::Real, MZ::MatrixZonotope)
    MZ.A0 .*= α
    @inbounds for i in ngens(MZ)
        MZ.Ai[i] .*= α
    end
    return MZ
end

"""
    *(a::Real, B::MatrixZonotope)

Alias to scale a matrix zonotope.
"""
@commutative Base.:*(a::Real, B::MatrixZonotope) = scale(a, B)
