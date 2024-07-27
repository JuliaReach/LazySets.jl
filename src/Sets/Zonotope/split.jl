function _split(Z::Zonotope, j::Int)
    c, G = Z.center, Z.generators

    c₁ = similar(c)
    G₁ = similar(G)
    Z₁ = Zonotope(c₁, G₁)

    c₂ = similar(c)
    G₂ = similar(G)
    Z₂ = Zonotope(c₂, G₂)

    return split!(Z₁, Z₂, Z, j)
end

function split!(Z₁::Zonotope, Z₂::Zonotope, Z::Zonotope, j::Int)
    c, G = Z.center, Z.generators
    n, p = size(G)
    @assert 1 <= j <= p "cannot split a zonotope with $p generators along " *
                        "index $j"

    c₁, G₁ = Z₁.center, Z₁.generators
    c₂, G₂ = Z₂.center, Z₂.generators
    copyto!(G₁, G)

    @inbounds for i in 1:n
        α = G[i, j] / 2
        c₁[i] = c[i] - α
        c₂[i] = c[i] + α
        G₁[i, j] = α
    end
    copyto!(G₂, G₁)

    return _split_ret(Z₁, Z₂)
end

_split_ret(Z₁::Zonotope, Z₂::Zonotope) = (Z₁, Z₂)

function load_StaticArraysCore_split()
    return quote
        using .StaticArraysCore: MMatrix, MVector

        function _split_ret(Z₁::Zonotope{N,SV,SM},
                            Z₂::Zonotope{N,SV,SM}) where {N,n,p,SV<:MVector{n,N},SM<:MMatrix{n,p,N}}
            Z₁ = Zonotope(StaticArraysCore.SVector{n}(Z₁.center),
                          StaticArraysCore.SMatrix{n,p}(Z₁.generators))
            Z₂ = Zonotope(StaticArraysCore.SVector{n}(Z₂.center),
                          StaticArraysCore.SMatrix{n,p}(Z₂.generators))
            return Z₁, Z₂
        end
    end
end  # load_StaticArraysCore_split

function _split(Z::Zonotope, gens::AbstractVector, n::AbstractVector)
    p = length(gens)
    @assert p == length(n) "the number of generators $(length(n)) does not " *
                           "match the number of indicated partitions $p"

    @assert p <= ngens(Z) "the number of generators to split ($p) is greater " *
                          "than the number of generators of the zonotope ($(ngens(Z)))"

    Zs = [Z]
    for (i, g) in enumerate(gens)
        for j in 1:n[i]
            km = length(Zs)
            for k in 1:km
                append!(Zs, split(Zs[k], g))
            end
            deleteat!(Zs, 1:km)
        end
    end
    return Zs
end
