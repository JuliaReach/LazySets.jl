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

function symmetric_interval_hull(H::Hyperrectangle{N}) where {N}
    n = dim(H)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        r[i] = _maxabs(H.center[i], H.radius[i])
    end
    return Hyperrectangle(zeros(N, n), r)
end

function symmetric_interval_hull(S::AbstractSingleton{N}) where {N}
    n = dim(S)
    r = abs.(element(S))
    return Hyperrectangle(zeros(N, n), r)
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
