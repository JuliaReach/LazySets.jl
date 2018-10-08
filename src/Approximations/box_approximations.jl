# ===================================
# Approximations in the infinity norm
# ===================================

"""
    box_approximation(S::LazySet; upper_bound::Bool=false)::Hyperrectangle

Overapproximate a convex set by a tight hyperrectangle.

### Input

- `S`           -- convex set
- `upper_bound` -- (optional, default: `false`) use overapproximation in support
                   function computation?

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet;
                           upper_bound::Bool=false
                          )::Hyperrectangle
    (c, r) = box_approximation_helper(S; upper_bound=upper_bound)
    return Hyperrectangle(c, r)
end

# special case: Hyperrectangle
box_approximation(S::Hyperrectangle; upper_bound::Bool=false) = S

# special case: other rectangle
box_approximation(S::AbstractHyperrectangle; upper_bound::Bool=false) =
    Hyperrectangle(center(S), radius_hyperrectangle(S))

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    box_approximation_symmetric(S::LazySet{N}; upper_bound::Bool=false
                               )::Hyperrectangle{N} where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle centered in the origin.

### Input

- `S`           -- convex set
- `upper_bound` -- (optional, default: `false`) use overapproximation in support
                   function computation?

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N};
                                     upper_bound::Bool=false
                                    )::Hyperrectangle{N} where {N<:Real}
    (c, r) = box_approximation_helper(S; upper_bound=upper_bound)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    symmetric_interval_hull

Alias for `box_approximation_symmetric`.
"""
symmetric_interval_hull = box_approximation_symmetric

"""
    box_approximation_helper(S::LazySet{N}; upper_bound::Bool=false
                            ) where {N<:Real}

Common code of `box_approximation` and `box_approximation_symmetric`.

### Input

- `S`           -- convex set
- `upper_bound` -- (optional, default: `false`) use overapproximation in support
                   function computation?

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
@inline function box_approximation_helper(S::LazySet{N};
                                          upper_bound::Bool=false
                                         ) where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    ρ_rec = upper_bound ? ρ_upper_bound : ρ
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ_rec(d, S)
        d[i] = -one_N
        hbottom = -ρ_rec(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
    end
    return c, r
end

"""
    ballinf_approximation(S::LazySet{N}; upper_bound::Bool=false
                         )::BallInf{N} where {N<:Real}

Overapproximate a convex set by a tight ball in the infinity norm.

### Input

- `S`           -- convex set
- `upper_bound` -- (optional, default: `false`) use overapproximation in support
                   function computation?

### Output

A tight ball in the infinity norm.

### Algorithm

The center and radius of the box are obtained by evaluating the support function
of the given convex set along the canonical directions.
"""
function ballinf_approximation(S::LazySet{N};
                               upper_bound::Bool=false
                              )::BallInf{N} where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    ρ_rec = upper_bound ? ρ_upper_bound : ρ
    n = dim(S)
    c = Vector{N}(undef, n)
    r = zero_N
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ_rec(d, S)
        d[i] = -one_N
        hbottom = -ρ_rec(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        rcur = (htop - hbottom) / 2
        if (rcur > r)
            r = rcur
        end
    end
    return BallInf(c, r)
end
