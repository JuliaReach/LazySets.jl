# ================================
# Concrete symmetric interval hull
# ================================

"""
    symmetric_interval_hull(S::LazySet{N}) where {N}

Overapproximate a set by a tight hyperrectangle centered in the origin.

### Input

- `S` -- set

### Output

A tight hyperrectangle that is centrally symmetric wrt. the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated in the canonical directions.

### Notes

The result is a hyperrectangle and hence in particular convex.

An alias for this function is `box_approximation_symmetric`.
"""
function symmetric_interval_hull(S::LazySet{N}) where {N}
    # fallback returns a hyperrectangular set
    (c, r) = box_approximation_helper(S)
    if r[1] < 0
        return EmptySet{N}(dim(S))
    end
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    box_approximation_symmetric

Alias for `symmetric_interval_hull`.
"""
box_approximation_symmetric = symmetric_interval_hull

# ===============
# Specializations
# ===============

# empty set specialization
symmetric_interval_hull(∅::EmptySet) = ∅

function LazySets.concretize(sih::SymmetricIntervalHull)
    return symmetric_interval_hull(LazySets.concretize(sih.X))
end

# interval specialization
function symmetric_interval_hull(x::Interval)
    abs_inf = abs(min(x))
    abs_sup = abs(max(x))
    bound = max(abs_sup, abs_inf)
    return Interval(-bound, bound)
end

# hyperrectangle specialization
@inline function _maxabs(c::N, r::N) where {N}
    return c >= zero(N) ? c + r : -c + r
end

function symmetric_interval_hull(H::AbstractHyperrectangle{N}) where {N}
    n = dim(H)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        r[i] = _maxabs(center(H, i), radius_hyperrectangle(H, i))
    end
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(S::AbstractSingleton{N}) where {N}
    n = dim(S)
    r = abs.(element(S))
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(P::Union{VPolygon{N}, VPolytope{N}}) where {N}
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

function symmetric_interval_hull(X::LinearMap{N, <:AbstractSingleton}) where {N}
    n = dim(X)
    r = abs.(X.M * element(X.X))
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(X::MinkowskiSum{N, <:AbstractSingleton, <:AbstractSingleton}) where {N}
    n = dim(X)
    r = abs.(element(X.X) + element(X.Y))
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(lm::LinearMap{N, <:AbstractHyperrectangle}) where {N}
    n = dim(lm)
    M = lm.M
    H = lm.X
    r = abs.(M * center(H)) .+ abs.(M) * radius_hyperrectangle(H)
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(E::ExponentialMap{N, <:AbstractSingleton}) where {N}
    require(:Expokit; fun_name="symmetric_interval_hull")

    v = expmv(one(N), E.spmexp.M, element(E.X))
    c = zeros(N, dim(E))
    r = abs.(v)
    return Hyperrectangle(c, r)
end

function symmetric_interval_hull(E::ExponentialMap{N, <:AbstractHyperrectangle}) where {N}
    require(:Expokit; fun_name="symmetric_interval_hull")

    H = set(E)
    n = dim(H)
    x = zeros(N, n)
    r = abs.(expmv(one(N), E.spmexp.M, center(H)))
    @inbounds for i in 1:n
        x[i] = radius_hyperrectangle(H, i)
        if isapproxzero(x[i])
            x[i] = 0
            continue
        end
        r .+= abs.(expmv(one(N), E.spmexp.M, x))
        x[i] = 0
    end
    return Hyperrectangle(zeros(N, n), r)
end
