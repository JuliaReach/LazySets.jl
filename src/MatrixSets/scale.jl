function scale(α::Real ,Z::MatrixZonotope)
    A0 = α * center(Z)
    Ai = [α* gen for gen in generators(Z)]
    return MatrixZonotope(A0, Ai)
end

"""
# Extended help

    scale!(α::Real, Z::Zonotope)

### Algorithm

In place scaling of a matrix zonotope `Z`.
"""
function scale!(α::Real, Z::MatrixZonotope)
    Z.A0 .*= α
    @inbounds for i in ngens(Z)
        Z.Ai[i] .*= α
    end
    return Z
end