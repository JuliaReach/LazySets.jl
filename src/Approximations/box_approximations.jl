# ===================================
# Approximations in the infinity norm
# ===================================

"""
    box_approximation(S::LazySet)::Hyperrectangle

Overapproximate a convex set by a tight hyperrectangle.

### Input

- `S`           -- convex set

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet{N};
                          )::Union{Hyperrectangle{N}, EmptySet{N}} where N<:Real
    (c, r) = box_approximation_helper(S)
    if r[1] < 0
        return EmptySet{N}()
    end
    return Hyperrectangle(c, r)
end

# special case: Hyperrectangle
box_approximation(S::Hyperrectangle) = S

# special case: other rectangle
box_approximation(S::AbstractHyperrectangle) =
    Hyperrectangle(center(S), radius_hyperrectangle(S))

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    box_approximation_symmetric(S::LazySet{N}
                               )::Union{Hyperrectangle{N}, EmptySet{N}}
                                where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle centered in the origin.

### Input

- `S`           -- convex set

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N};
                                    )::Union{Hyperrectangle{N},
                                             EmptySet{N}} where {N<:Real}
    (c, r) = box_approximation_helper(S)
    if r[1] < 0
        return EmptySet{N}()
    end
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    symmetric_interval_hull

Alias for `box_approximation_symmetric`.
"""
symmetric_interval_hull = box_approximation_symmetric

"""
    box_approximation_helper(S::LazySet{N};
                            ) where {N<:Real}

Common code of `box_approximation` and `box_approximation_symmetric`.

### Input

- `S`           -- convex set

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
                                         ) where {N<:Real}
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
    ballinf_approximation(S::LazySet{N};
                         )::BallInf{N} where {N<:Real}

Overapproximate a convex set by a tight ball in the infinity norm.

### Input

- `S`           -- convex set

### Output

A tight ball in the infinity norm.

### Algorithm

The center and radius of the box are obtained by evaluating the support function
of the given convex set along the canonical directions.
"""
function ballinf_approximation(S::LazySet{N};
                              )::Union{BallInf{N}, EmptySet{N}} where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    n = dim(S)
    c = Vector{N}(undef, n)
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
        elseif rcur < 0
            # contradicting bounds => set is empty
            return EmptySet{N}()
        end
    end
    return BallInf(c, r)
end
