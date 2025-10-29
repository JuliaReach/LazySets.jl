"""
    quadratic_map(Q::Vector{<:AbstractMatrix}, S::SimpleSparsePolynomialZonotope)

Return the quadratic map of a simple sparse polynomial zonotope.

### Input

- `Q` -- vector of square matrices
- `S` -- simple sparse polynomial zonotope

### Output

The quadratic map of `P` represented as a simple sparse polynomial zonotope.

### Algorithm

This method implements [KochdumperA21; Proposition 12](@citet).
See also [Kochdumper21a; Proposition 3.1.30](@citet).
"""
function quadratic_map(Q::Vector{<:AbstractMatrix}, S::SimpleSparsePolynomialZonotope)
    m = length(Q)
    c = center(S)
    h = ngens(S)
    G = genmat(S)
    E = expmat(S)

    cnew = similar(c, m)
    Gnew = similar(G, m, h^2 + h)
    QiG = similar(Q)
    @inbounds for (i, Qi) in enumerate(Q)
        cnew[i] = dot(c, Qi, c)
        Gnew[i, 1:h] = c' * (Qi + Qi') * G
        QiG[i] = Qi * G
    end

    Enew = repeat(E, 1, h + 1)
    @inbounds for i in 1:h
        idxstart = h * i + 1
        idxend = (i + 1) * h
        Enew[:, idxstart:idxend] .+= E[:, i]
        for j in eachindex(QiG)
            Gnew[j, idxstart:idxend] = G[:, i]' * QiG[j]
        end
    end
    Z = SimpleSparsePolynomialZonotope(cnew, Gnew, Enew)
    return remove_redundant_generators(Z)
end

"""
    quadratic_map(Q::Vector{<:AbstractMatrix}, S1::SimpleSparsePolynomialZonotope,
                  S2::SimpleSparsePolynomialZonotope)

Return the quadratic map of two simple sparse polynomial zonotopes.
The quadratic map is the set
```math
    \\{x \\mid xᵢ = s₁ᵀQᵢs₂, s₁ ∈ S₁, s₂ ∈ S₂, Qᵢ ∈ Q\\}.
```

### Input

- `Q`  -- vector of square matrices
- `S1` -- simple sparse polynomial zonotope
- `S2` -- simple sparse polynomial zonotope

### Output

The quadratic map of the given simple sparse polynomial zonotopes represented as
a simple sparse polynomial zonotope.

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.30](@citet).
"""
function quadratic_map(Q::Vector{<:AbstractMatrix}, S1::SimpleSparsePolynomialZonotope,
                       S2::SimpleSparsePolynomialZonotope)
    @assert nparams(S1) == nparams(S2) "the number of parameters must be equal"

    c1 = center(S1)
    c2 = center(S2)
    G1 = genmat(S1)
    G2 = genmat(S2)
    E1 = expmat(S1)
    E2 = expmat(S2)

    c = [dot(c1, Qi, c2) for Qi in Q]

    Ghat1 = reduce(vcat, c2' * Qi' * G1 for Qi in Q)
    Ghat2 = reduce(vcat, c1' * Qi * G2 for Qi in Q)

    Gbar = reduce(hcat, reduce(vcat, gj' * Qi * G2 for Qi in Q) for gj in eachcol(G1))
    Ebar = reduce(hcat, E2 .+ e1j for e1j in eachcol(E1))

    G = hcat(Ghat1, Ghat2, Gbar)
    E = hcat(E1, E2, Ebar)

    return remove_redundant_generators(SimpleSparsePolynomialZonotope(c, G, E))
end
