# ===========================================================
# Concrete overapproximation with a ball in the infinity norm
# ===========================================================

"""
    ballinf_approximation(S::LazySet)

Overapproximate a set by a tight ball in the infinity norm.

### Input

- `S` -- set

### Output

A tight ball in the infinity norm.

### Algorithm

The center and radius of the box are obtained by evaluating the support function
of the given set along the canonical directions.
"""
function ballinf_approximation(S::LazySet{N}) where {N}
    n = dim(S)
    c = Vector{N}(undef, n)
    r = zero(N)
    d = zeros(N, n)

    @inbounds for i in 1:n
        d[i] = one(N)
        htop = ρ(d, S)
        d[i] = -one(N)
        hbottom = -ρ(d, S)
        d[i] = zero(N)
        c[i] = (htop + hbottom) / 2
        rcur = (htop - hbottom) / 2
        if (rcur > r)
            r = rcur
        elseif !_geq(rcur, zero(N))
            # contradicting bounds => set is empty
            return EmptySet{N}(dim(S))
        end
    end
    return BallInf(c, r)
end

# ===============
# Specializations
# ===============

# empty set specialization
ballinf_approximation(∅::EmptySet) = ∅
