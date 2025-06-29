function rectify(U::Universe)
    require(@__MODULE__, :LazySets; fun_name="linear_map")

    N = eltype(U)
    n = dim(U)
    clist = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, n)
    @inbounds for i in 1:n
        clist[i] = HalfSpace(SingleEntryVector(i, n, N(-1)), zero(N))
    end
    return HPolyhedron(clist)
end
