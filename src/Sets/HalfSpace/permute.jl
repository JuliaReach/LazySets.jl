@validate function permute(H::HalfSpace, p::AbstractVector{Int})
    return HalfSpace(H.a[p], H.b)
end
