import Base.LinAlg:norm
import LazySets:radius,  # visible outside module Approximations
                diameter

# ===================================
# Approximations in the infinity norm
# ===================================

"""
    box_approximation(S::LazySet)::Hyperrectangle

Overapproximate a convex set by a tight hyperrectangle.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet)::Hyperrectangle
    (c, r) = box_approximation_helper(S)
    return Hyperrectangle(c, r)
end

box_approximation(S::Hyperrectangle) = S
box_approximation(S::BallInf) = Hyperrectangle(S.center, fill(S.radius, dim(S)))

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    box_approximation_symmetric(S::LazySet{N})::Hyperrectangle{N} where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle centered in the origin.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N}
                                    )::Hyperrectangle{N} where {N<:Real}
    (c, r) = box_approximation_helper(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    symmetric_interval_hull

Alias for `box_approximation_symmetric`.
"""
symmetric_interval_hull = box_approximation_symmetric

"""
    box_approximation_helper(S::LazySet)

Common code of `box_approximation` and `box_approximation_symmetric`.

### Input

- `S` -- convex set

### Output

A tuple containing the data that is needed to construct a tightly
overapproximating hyperrectangle.

- `c` -- center
- `r` -- radius

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given convex set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.
"""
@inline function box_approximation_helper(S::LazySet{N}) where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    n = dim(S)
    c = Vector{N}(n)
    r = Vector{N}(n)
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ(d, S)
        d[i] = -one_N
        hbottom = -ρ(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
    end
    return c, r
end

"""
    ballinf_approximation(S::LazySet{N})::BallInf{N} where {N<:Real}

Overapproximate a convex set by a tight ball in the infinity norm.

### Input

- `S` -- convex set

### Output

A tight ball in the infinity norm.

### Algorithm

The center and radius of the box are obtained by evaluating the support function
of the given convex set along the canonical directions.
"""
function ballinf_approximation(S::LazySet{N})::BallInf{N} where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    n = dim(S)
    c = Vector{N}(n)
    r = zero_N
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ(d, S)
        d[i] = -one_N
        hbottom = -ρ(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        rcur = (htop - hbottom) / 2
        if (rcur > r)
            r = rcur
        end
    end
    return BallInf(c, r)
end

# =======================================================
# Metric properties of sets computed using Approximations
# =======================================================

"""
    norm(S::LazySet, [p]::Real=Inf)

Return the norm of a convex set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `S` -- convex set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(S::LazySet, p::Real=Inf)
    if p == Inf
        return norm(ballinf_approximation(S), p)
    else
        error("the norm for this value of p=$p is not implemented")
    end
end

"""
    radius(S::LazySet, [p]::Real=Inf)

Return the radius of a convex set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `S` -- convex set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.
"""
function radius(S::LazySet, p::Real=Inf)
    if p == Inf
        return radius(ballinf_approximation(S)::BallInf, p)
    else
        error("the radius for this value of p=$p is not implemented")
    end
end

"""
    diameter(S::LazySet, [p]::Real=Inf)

Return the diameter of a convex set.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `S` -- convex set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
function diameter(S::LazySet, p::Real=Inf)
    if p == Inf
        return radius(S, p) * 2
    else
        error("the diameter for this value of p=$p is not implemented")
    end
end
