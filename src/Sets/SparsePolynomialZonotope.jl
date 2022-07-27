export SparsePolynomialZonotope, expmat, nparams,
       linear_map, quadratic_map, remove_redundant_generators

struct SparsePolynomialZonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}, VI<:AbstractVector{<:Integer}} <: AbstractPolynomialZonotope{N}
    c::VN
    G::MN
    GI::MN
    E::ME
    idx::VI

    function SparsePolynomialZonotope(c::VN, G::MN, GI::MN, E::ME, idx::VI) where {N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ME<:AbstractMatrix{<:Integer}}
        @assert length(c) == size(G, 1) throw(DimensionMismatch("c and G should have the same number of rows"))
        @assert length(c) == size(GI, 1) throw(DimensionMismatch("c and GI should have the same number of rows"))
        @assert size(G, 2) == size(E, 2) throw(DimensionMismatch("G and E should have the same number of columns"))
        @assert all(>=(0), E) throw(ArgumentError("E should have non-negative integers"))

        return new{N, VN, MN, ME}(c, G, GI, E, VI)
    end
end


const SPZ = SimpleSparsePolynomialZonotope

dim(P::SPZ) = size(P.c, 1)

ndependentgens(P::SPZ) = size(P.G, 2)

nindependentgens(P::SPZ) = size(P.GI, 2)

ngens(P::SPZ) = ndependentgens(P) + nindependentgens(P)

nparams(P::SPZ) = size(P.E, 1)

order(P::SPZ) = ngens(P) // dim(P)

center(P::SPZ) = P.c

dependent_genmat(P::SPZ) = P.G

independent_genmat(P::SPZ) = P.GI

expmat(P::SPZ) = P.E

indexvector(P::SPZ) = P.idx


function linear_map(M::AbstractMatrix, P::SPZ)
    return SparsePolynomialZonotope(M * center(P),
                                    M * dependent_genmat(P),
                                    M * independent_genmat(P),
                                    expmat(P),
                                    indexvector(P)
                                    )
end

function quadratic_map(Q::Vector{MT}, S::SparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
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
    return SparsePolynomialZonotope(cnew, Gnew, Enew)
end

function remove_redundant_generators(S::SparsePolynomialZonotope)

    c = center(S)
    G = genmat(S)
    E = expmat(S)

    Gnew = Matrix{eltype(G)}(undef, size(G, 1), 0)
    Enew = Matrix{eltype(E)}(undef, size(E, 1), 0)
    cnew = copy(c)

    visited_exps = Dict{Vector{Int}, Int}()
    @inbounds for (gi, ei) in zip(eachcol(G), eachcol(E))
        iszero(gi) && continue
        if iszero(ei)
            cnew += gi
        elseif haskey(visited_exps, ei) # repeated exponent
            idx = visited_exps[ei]
            Gnew[:, idx] += gi
        else
            Gnew = hcat(Gnew, gi)
            Enew = hcat(Enew, ei)
            visited_exps[ei] = size(Enew, 2)
        end
    end

    return SparsePolynomialZonotope(cnew, Gnew, Enew)
end
