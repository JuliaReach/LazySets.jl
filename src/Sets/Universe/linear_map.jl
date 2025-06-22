function linear_map(M::AbstractMatrix, U::Universe)
    @assert size(M, 2) == dim(U) "cannot apply a $(size(M))-dimensional " *
                                 "matrix to a $(dim(U))-dimensional set"

    N = promote_type(eltype(M), eltype(U))

    # a zero row in M is a projection to zero in that dimension
    clist = Vector{HalfSpace{N,SingleEntryVector{N}}}()
    m = size(M, 1)
    for (j, cj) in enumerate(eachrow(M))
        if iszero(cj)
            push!(clist, HalfSpace(SingleEntryVector(j, m, N(1)), N(0)))
            push!(clist, HalfSpace(SingleEntryVector(j, m, N(-1)), N(0)))
        end
    end
    if !isempty(clist)
        return HPolyhedron(clist)
    end

    return Universe{N}(m)
end
