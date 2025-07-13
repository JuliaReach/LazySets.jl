@validate function permute(S::Singleton, p::AbstractVector{Int})
    e = S.element[p]
    return Singleton(e)
end
