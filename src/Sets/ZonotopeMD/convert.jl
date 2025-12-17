function convert(::Type{ZonotopeMD}, Z::AbstractZonotope)
    N = eltype(Z)
    n = dim(Z)
    Gvec = Vector{N}[]
    d = zeros(N, n)
    for g in generators(Z)
        i = find_unique_nonzero_entry(g)
        if i > 0
            d[i] += abs(g[i])
        else
            push!(Gvec, g)
        end
    end
    G = reduce(hcat, Gvec; init=zeros(N, n, 0))
    return ZonotopeMD(center(Z), G, d)
end
