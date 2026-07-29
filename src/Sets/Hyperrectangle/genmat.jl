function genmat(H::Hyperrectangle{N,<:AbstractVector,<:SparseVector{N}}) where {N}
    n = dim(H)
    nze, nzv = findnz(H.radius)
    return sparse(nze, 1:length(nze), nzv, n, length(nze))
end
