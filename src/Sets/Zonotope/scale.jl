function scale(α::Real, Z::Zonotope)
    if iszero(α)
        # generator matrix would consist of only zero columns -> create new matrix with no columns
        N = promote_type(typeof(α), eltype(Z))
        return Zonotope(zeros(N, dim(Z)), zeros(N, dim(Z), 0))
    else
        return Zonotope(α * Z.center, α * Z.generators)
    end
end

"""
# Extended help

    scale!(α::Real, Z::Zonotope)

### Algorithm

The result is obtained by applying the numerical scale to the center and
generators.
"""
function scale!(α::Real, Z::Zonotope)
    Z.center .*= α
    Z.generators .*= α
    return Z
end
