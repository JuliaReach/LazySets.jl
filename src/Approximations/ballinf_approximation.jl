# ===========================================================
# Concrete overapproximation with a ball in the infinity norm
# ===========================================================

"""
    ballinf_approximation(S::LazySet)

Overapproximate a convex set by a tight ball in the infinity norm.

### Input

- `S` -- convex set

### Output

A tight ball in the infinity norm.
"""
function ballinf_approximation(S::LazySet)
    return overapproximate(S, BallInf)
end

# ===============
# Specializations
# ===============

# empty set specialization
ballinf_approximation(∅::EmptySet) = ∅
