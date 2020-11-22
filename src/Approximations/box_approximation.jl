# ================================================
# Concrete overapproximation with a hyperrectangle
# ================================================

"""
    box_approximation(S::LazySet{N}) where {N<:Real}

Overapproximate a set by a tight hyperrectangle.

### Input

- `S` -- set

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
    box_approximation_helper(S::LazySet{N}) where {N<:Real}

Common code of `box_approximation` and `box_approximation_symmetric`.

### Input

- `S` -- set

### Output

A tuple containing the data that is needed to construct a tightly
overapproximating hyperrectangle.

- `c` -- center
- `r` -- radius

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.
"""
@inline function box_approximation_helper(S::LazySet{N}) where {N<:Real}
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        htop = ρ(SingleEntryVector(i, n, one(N)), S)
        hbottom = -ρ(SingleEntryVector(i, n, -one(N)), S)
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
