function extrema(P::SparsePolynomialZonotope; algorithm::String="zonotope")
    if algorithm == "zonotope"
        return _extrema_polyzono_zonotope(P)
    elseif algorithm == "lowhigh"
        return _extrema_lowhigh(P)
    else
        throw(ArgumentError("unknown algorithm \"$algorithm\""))
    end
end

# See [KochdumperSAB23; Proposition 1](@citet).
function _extrema_polyzono_zonotope(P::SparsePolynomialZonotope{N}) where {N}
    G = genmat_dep(P)
    E = expmat(P)
    n = dim(P)
    g₁ = zeros(N, n)
    g₂ = zeros(N, n)
    g₃ = zeros(N, n)
    for (j, Gj) in enumerate(eachcol(G))
        if all(iseven.(@view E[:, j]))
            g₁ .+= Gj
            g₂ .+= abs.(Gj)
        else
            g₃ .+= abs.(Gj)
        end
    end
    g₄ = zeros(N, n)
    for GIj in eachcol(genmat_indep(P))
        g₄ .+= abs.(GIj)
    end
    c = center(P) + N(1 / 2) .* g₁
    r = N(1 / 2) .* g₂ + g₃ + g₄
    return (c .- r, c .+ r)
end
