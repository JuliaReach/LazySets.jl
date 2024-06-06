function complement(X::Interval)
    N = eltype(X)
    L = HalfSpace(SingleEntryVector(1, 1, one(N)), min(X))
    H = HalfSpace(SingleEntryVector(1, 1, -one(N)), -max(X))
    return UnionSet(L, H)
end
