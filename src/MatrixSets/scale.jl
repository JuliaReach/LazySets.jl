function scale(α::Real, Z::MatrixZonotope)
    return scale!(α, copy(Z))
end

function scale!(α::Real, Z::MatrixZonotope)
    Z.A0 .*= α
    @inbounds for i in ngens(Z)
        Z.Ai[i] .*= α
    end
    return Z
end
