"""
See Def. 4 in [1].

[1] Luo, Ertai, Niklas Kochdumper, and Stanley Bak. "Reachability analysis for linear systems with uncertain parameters using polynomial zonotopes." Proceedings of the 26th ACM International Conference on Hybrid Systems: Computation and Control. 2023.
"""
struct MatrixZonotope{N,MN<:AbstractMatrix{N}}
    "Center matrix."
    A0::MN
    "Generator matrices."
    Ai::Vector{MN}
    "Unique identifiers."
    id::Vector{Int64}
    function MatrixZonotope(A0::MN, Ai::Vector{MN},
                            id::Vector{Int64}=collect(1:length(Ai))) where {N,MN<:AbstractMatrix{N}}
        length(Ai) == length(unique(id)) ||
            throw(ArgumentError("The number of generator matrices doesn't match the id's."))
        if length(Ai) > 0
            (allequal(size.(Ai)) && size(A0) == first(Ai)) ||
                throw(ArgumentError("The size of all generator matrices should match."))
        end
        return new{N,MN}(A0, Ai, id)
    end
end

"""
See Prop. 1 in [1].

[1] Luo, Ertai, Niklas Kochdumper, and Stanley Bak. "Reachability analysis for linear systems with uncertain parameters using polynomial zonotopes." Proceedings of the 26th ACM International Conference on Hybrid Systems: Computation and Control. 2023.
"""
function linear_map(MZ::MatrixZonotope, P::SimpleSparsePolynomialZonotope)
    A0 = MZ.A0
    Ai = MZ.Ai
    @assert size(A0, 2) == dim(P)
    cZ = center(P)
    GZ = genmat(P)
    EZ = expmat(P)
    w = length(MZ.Ai)

    c = A0 * cZ
    G = vcat(A0 * GZ, [A * cZ for A in Ai], [A * GZ for A in Ai])
    Iw = Diagonal(ones(Int, w))
    E₁, E₂, id = merge_id(P.id, MZ.id, EZ, Iw)
    h = size(GZ, 2)
    Bh = ones(1, h)
    E = vcat(E₁, E₂, [E₂[:, j] * Bh + E₁ for j in 1:w])
    # TODO Make SimpleSparsePolynomialZonotope have an `id` field.
    return SimpleSparsePolynomialZonotope(c, G, E)
end

function merge_id(id1::Vector{Int64}, id2::Vector{Int64}, E₁::AbstractMatrix{N},
                  E₂::AbstractMatrix{N}) where {N}
    p₁, h₁ = size(E₁)
    p₂, h₂ = size(E₂)
    if !(length(id1) == p₁ && length(id2) == p₂)
        throw(ArgumentError("Incompatible dimensions."))
    end

    # Indices of elements of id2 which do not belong to id1.
    K = findall(∈(setdiff(id2, id1)), id2)
    k = length(K)

    id = vcat(id1, id2[K])
    Ē₁ = vcat(E₁, zeros(N, k, h₁))

    Ē₂ = zeros(N, p₁ + k, h₂)
    for i in 1:(p₁ + k)
        j = findfirst(==(id[i]), id₂)
        if !isnothing(j)
            Ē₂[i, :] = E₂[j, :]
        end
    end
    return Ē₁, Ē₂, id
end
