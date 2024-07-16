function complement(X::Interval)
    require(@__MODULE__, :LazySets; fun_name="complement")

    N = eltype(X)
    L = HalfSpace(SingleEntryVector(1, 1, one(N)), min(X))
    H = HalfSpace(SingleEntryVector(1, 1, -one(N)), -max(X))
    return UnionSet(L, H)
end
