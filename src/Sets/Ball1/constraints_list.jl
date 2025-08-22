"""
# Extended help

    constraints_list(P::Ball1)

### Notes

In ``n`` dimensions there are ``2^n`` constraints (unless the radius is 0).

### Algorithm

The constraints can be defined as ``d_i^T (x-c) ≤ r`` for all ``d_i``, where
``d_i`` is a vector with elements ``1`` or ``-1`` in ``n`` dimensions. To span
all possible ``d_i``, the function `Iterators.product` is used.
"""
@validate function constraints_list(B::Ball1)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    n = dim(B)
    c, r = B.center, B.radius
    N = eltype(B)
    clist = Vector{HalfSpace{N,Vector{N}}}(undef, 2^n)
    for (i, di) in enumerate(Iterators.product([[one(N), -one(N)] for i in 1:n]...))
        d = collect(di) # tuple -> vector
        clist[i] = HalfSpace(d, dot(d, c) + r)
    end
    return clist
end
