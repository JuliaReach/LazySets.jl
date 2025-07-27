function scale(α::Real, MZ::MatrixZonotope)
    return scale!(α, copy(MZ))
end

function scale!(α::Real, MZ::MatrixZonotope)
    MZ.A0 .*= α
    @inbounds for i in 1:ngens(MZ)
        MZ.Ai[i] .*= α
    end
    return MZ
end

"""
    *(a::Real, B::MatrixZonotope)

Alias to scale a matrix zonotope.
"""
@commutative Base.:*(α::Real, MZ::MatrixZonotope) = scale(α, MZ)
