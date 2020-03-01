# ===================================
# Approximations in the infinity norm
# ===================================

"""
    box_approximation(S::LazySet{N}) where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle.
"""
function box_approximation(S::LazySet{N}) where {N<:Real}
    return overapproximate(S, Hyperrectangle)
end

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    box_approximation_symmetric(S::LazySet{N}) where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle centered in the origin.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N}) where {N<:Real}
    # fallback returns a hyperrectangular set
    (c, r) = box_approximation_helper(S)
    if r[1] < 0
        return EmptySet{N}(dim(S))
    end
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

# special case: empty set
box_approximation_symmetric(∅::EmptySet) = ∅

"""
    symmetric_interval_hull

Alias for `box_approximation_symmetric`.
"""
symmetric_interval_hull = box_approximation_symmetric

# interval specialization
function box_approximation_symmetric(x::Interval)
    abs_inf = abs(min(x))
    abs_sup = abs(max(x))
    bound = max(abs_sup, abs_inf)
    return Interval(-bound, bound)
end

# hyperrectangle specialization
@inline function _maxabs(c::N, r::N) where {N}
    return c >= zero(N) ? c + r : -c + r
end

function box_approximation_symmetric(H::Hyperrectangle{N}) where {N<:Real}
    n = dim(H)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        r[i] = _maxabs(H.center[i], H.radius[i])
    end
    return Hyperrectangle(zeros(N, n), r)
end

"""
    box_approximation_helper(S::LazySet{N}) where {N<:Real}

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
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ(d, S)
        d[i] = -one_N
        hbottom = -ρ(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
        if r[i] < 0
            # contradicting bounds => set is empty
            # terminate with first radius entry being negative
            r[1] = r[i]
            break
        end
    end
    return c, r
end

"""
    ballinf_approximation(S::LazySet{N}) where {N<:Real}

Overapproximate a convex set by a tight ball in the infinity norm.

### Input

- `S` -- convex set

### Output

A tight ball in the infinity norm.
"""
function ballinf_approximation(S::LazySet{N}) where {N<:Real}
    return overapproximate(S, BallInf)
end

# special case: empty set
ballinf_approximation(∅::EmptySet) = ∅
