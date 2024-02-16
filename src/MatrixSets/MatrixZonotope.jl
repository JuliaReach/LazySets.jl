struct MatrixZonotope{N,MN<:AbstractMatrix{N}}
    A0::MN
    Ai::Vector{MN}
end

function linear_map(MZ::MatrixZonotope, P::SimpleSparsePolynomialZonotope)
    @assert size(MZ.A0, 2) == dim(P)
    cZ = center(P)
    GZ = genmat(P)
    EZ = expmat(P)
    w = length(MZ.Ai)

    c = MZ.A0 * cZ
    G = vcat(MZ.A0 * GZ, [A * cZ for A in MZ.Ai], [A * GZ for A in MZ.Ai])
    E₁, E₂ = merge_id(EZ, Diagonal(ones(Int, w)))
    E = vcat(E₁, E₂, [TODO...])
    return SimpleSparsePolynomialZonotope(c, G, E)
end

# this implementation assumes that the IDs are 1:p₁ and 1:p₂
function merge_id(E₁, E₂)
    p₁ = size(E₁, 1)
    p₂ = size(E₂, 1)
    if p₁ < p₂
        E₁ = vcat(E₁, zeros(Int, p₂ - p₁))
    elseif p₁ > p₂
        E₂ = vcat(E₂, zeros(Int, p₁ - p₂))
    end
    return E₁, E₂
end
