"""
    symmetric_interval_hull(S::LazySet{N}) where {N}

Overapproximate a set by a tight hyperrectangle centered in the origin.

### Input

- `S` -- set

### Output

A tight hyperrectangle that is centrally symmetric wrt. the origin.

### Notes

The convenience alias `box_approximation_symmetric` is also available.

### Algorithm

The center of the box is the origin, and the radius is obtained via the
`extrema` function.
"""
function symmetric_interval_hull(S::LazySet{N}) where {N}
    n = dim(S)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        lo, hi = extrema(S, i)
        if !_geq(hi, lo)
            # contradicting bounds => set is empty
            return EmptySet{N}(n)
        end
        r[i] = max(abs(lo), abs(hi))
    end
    return Hyperrectangle(zeros(N, length(r)), r)
end

box_approximation_symmetric = symmetric_interval_hull

# concretization of a lazy SymmetricIntervalHull
function LazySets.concretize(sih::SymmetricIntervalHull)
    return symmetric_interval_hull(sih.X)
end

function symmetric_interval_hull(H::AbstractHyperrectangle{N}) where {N}
    n = dim(H)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        ci = center(H, i)
        ri = radius_hyperrectangle(H, i)
        r[i] = ci >= zero(N) ? ci + ri : -ci + ri
    end
    return Hyperrectangle(zeros(N, n), r)
end

symmetric_interval_hull(∅::EmptySet) = ∅

function symmetric_interval_hull(S::AbstractSingleton{N}) where {N}
    n = dim(S)
    r = abs.(element(S))
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(P::Union{VPolygon{N},VPolytope{N}}) where {N}
    n = dim(P)
    r = zeros(N, n)
    @inbounds for v in vertices(P)
        for i in 1:n
            r[i] = max(r[i], abs(v[i]))
        end
    end
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(L::LineSegment{N}) where {N}
    r = @inbounds [max(abs(L.p[1]), abs(L.q[1])),
                   max(abs(L.p[2]), abs(L.q[2]))]
    return Hyperrectangle(zeros(N, 2), r)
end

function symmetric_interval_hull(lm::LinearMap{N,<:AbstractSingleton}) where {N}
    n = dim(lm)
    r = abs.(lm.M * element(lm.X))
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(X::MinkowskiSum{N,<:AbstractSingleton,<:AbstractSingleton}) where {N}
    n = dim(X)
    r = abs.(element(first(X)) + element(second(X)))
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(lm::LinearMap{N,<:AbstractHyperrectangle}) where {N}
    n = dim(lm)
    M = lm.M
    H = lm.X
    r = abs.(M * center(H)) .+ abs.(M) * radius_hyperrectangle(H)
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(E::ExponentialMap{N,<:AbstractSingleton};
                                 backend=get_exponential_backend()) where {N}
    v = _expmv(backend, one(N), E.spmexp.M, element(E.X))
    c = zeros(N, dim(E))
    r = abs.(v)
    return Hyperrectangle(c, r)
end

function symmetric_interval_hull(E::ExponentialMap{N,<:AbstractHyperrectangle};
                                 backend=get_exponential_backend()) where {N}
    H = set(E)
    n = dim(H)
    x = zeros(N, n)
    r = abs.(_expmv(backend, one(N), E.spmexp.M, center(H)))
    @inbounds for i in 1:n
        x[i] = radius_hyperrectangle(H, i)
        if isapproxzero(x[i])
            x[i] = 0
            continue
        end
        r .+= abs.(_expmv(backend, one(N), E.spmexp.M, x))
        x[i] = 0
    end
    return Hyperrectangle(zeros(N, n), r)
end
