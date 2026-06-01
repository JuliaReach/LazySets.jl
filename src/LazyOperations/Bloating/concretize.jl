function concretize(B::Bloating{N}) where {N}
    return minkowski_sum(concretize(B.X), Ballp(B.p, zeros(N, dim(B)), B.ε))
end
