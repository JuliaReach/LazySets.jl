"""
    ρ(d::AbstractVector, L::Line)

Evaluate the support function of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, L::Line)
    if isapproxzero(dot(d, L.d))
        return dot(d, L.p)
    else
        N = eltype(L)
        return N(Inf)
    end
end
