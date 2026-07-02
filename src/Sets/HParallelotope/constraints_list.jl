function constraints_list(P::HParallelotope{N,VN}) where {N,VN}
    D, c = P.directions, P.offset
    return _constraints_list_hparallelotope(D, c, N, VN)
end

# see ext/LazySets/LazySetsHParallelotopeExt.jl
_constraints_list_hparallelotope(D, c, N, VN) = error()
