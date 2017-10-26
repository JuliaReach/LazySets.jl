# =======================================
# Support function of an abstract LazySet
# =======================================

"""
    ρ(d::Vector{Float64}, sf::LazySet)::Float64

Evaluate the support function of a set in a given direction.

### Input

- `d`  -- a real vector, the direction investigated
- `sf` -- a convex set

### Output

- `ρ(d, sf)` -- the support function
"""
function ρ(d::Vector{Float64}, sf::LazySet)::Float64
    return dot(d, σ(d, sf))
end

# alias
support_function = ρ
support_vector = σ
