"""
    overapproximate(S::LazySet{N},
                    ::Type{<:HPolygon},
                    [ɛ]::Real=Inf)::HPolygon where {N<:Real}

Return an approximation of a given 2D convex set.
If no error tolerance is given, or is `Inf`, the result is a box-shaped polygon.
Otherwise the result is an ɛ-close approximation as a polygon.

### Input

- `S` -- convex set, assumed to be two-dimensional
- `HPolygon` for dispatch
- `ɛ` -- (optional, default: `Inf`) error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(S::LazySet{N},
                         ::Type{<:HPolygon},
                         ɛ::Real=Inf)::HPolygon where {N<:Real}
    @assert dim(S) == 2
    if ɛ == Inf
        pe = σ(DIR_EAST(N), S)
        pn = σ(DIR_NORTH(N), S)
        pw = σ(DIR_WEST(N), S)
        ps = σ(DIR_SOUTH(N), S)
        constraints = Vector{LinearConstraint{eltype(pe)}}(4)
        constraints[1] = LinearConstraint(DIR_EAST(N), dot(pe, DIR_EAST(N)))
        constraints[2] = LinearConstraint(DIR_NORTH(N), dot(pn, DIR_NORTH(N)))
        constraints[3] = LinearConstraint(DIR_WEST(N), dot(pw, DIR_WEST(N)))
        constraints[4] = LinearConstraint(DIR_SOUTH(N), dot(ps, DIR_SOUTH(N)))
        return HPolygon(constraints)
    else
        return tohrep(approximate(S, ɛ))
    end
end

"""
    overapproximate(S::LazySet, ɛ::Real)::HPolygon

Alias for `overapproximate(S, HPolygon, ɛ)`.
"""
overapproximate(S::LazySet, ɛ::Real)::HPolygon = overapproximate(S, HPolygon, ɛ)

"""
    overapproximate(S::LazySet, Type{<:Hyperrectangle})::Hyperrectangle

Return an approximation of a given set as a hyperrectangle.

### Input

- `S` -- set
- `Hyperrectangle` for dispatch

### Output

A hyperrectangle.
"""
overapproximate(S::LazySet, ::Type{<:Hyperrectangle}) = box_approximation(S)

"""
    overapproximate(S::LazySet)::Hyperrectangle

Alias for `overapproximate(S, Hyperrectangle)`.
"""
overapproximate(S::LazySet)::Hyperrectangle = overapproximate(S, Hyperrectangle)

"""
    overapproximate(S::ConvexHull{N, Zonotope{N}, Zonotope{N}},
                     ::Type{<:Zonotope})::Zonotope where {N<:Real}

Overapproximate the convex hull of two zonotopes.

### Input

- `S` -- convex hull of two zonotopes of the same order
- `Zonotope` for dispatch

### Algorithm

This function implements the method proposed in
*Reachability of Uncertain Linear Systems Using Zonotopes, A. Girard, HSCC 2005*.
The convex hull of two zonotopes of the same order, that we write
``Z_j = ⟨c^{(j)}, g^{(j)}_1, …, g^{(j)}_p⟩`` for ``j = 1, 2``, can be
overapproximated as follows:

```math
CH(Z_1, Z_2) ⊆ \\frac{1}{2}⟨c^{(1)}+c^{(2)}, g^{(1)}_1+g^{(2)}_1, …, g^{(1)}_p+g^{(2)}_p, c^{(1)}-c^{(2)}, g^{(1)}_1-g^{(2)}_1, …, g^{(1)}_p-g^{(2)}_p⟩.
```

It should be noted that the output zonotope is not necessarily the minimal enclosing
zonotope, which is in general expensive in high dimensions. This is further investigated
in: *Zonotopes as bounding volumes, L. J. Guibas et al, Proc. of Symposium on Discrete Algorithms, pp. 803-812*.
"""
function overapproximate(S::ConvexHull{N, Zonotope{N}, Zonotope{N}},
                         ::Type{<:Zonotope})::Zonotope where {N<:Real}
    Z1, Z2 = S.X, S.Y
    @assert order(Z1) == order(Z2)
    center = (Z1.center+Z2.center)/2
    generators = [(Z1.generators .+ Z2.generators) (Z1.center - Z2.center) (Z1.generators .- Z2.generators)]/2
    return Zonotope(center, generators)
end

@require IntervalArithmetic begin

"""
    overapproximate(::LazySet{N}, ::Type{LazySets.Interval}) where {N<:Real}

Return the overapproximation of a real unidimensional set with an interval.

### Input

- `S` -- one-dimensional set
- `Interval` for dispatch

### Output

An interval.
"""
function overapproximate(S::LazySet{N}, ::Type{LazySets.Interval}) where {N<:Real}
    @assert dim(S) == 1
    lo = σ([-one(N)], S)[1]
    hi = σ([one(N)], S)[1]
    return LazySets.Interval(lo, hi)
end

end # @require

"""
    overapproximate(X::S, ::Type{S}) where {S <: LazySet}

Overapproximating a set of type `S` with type `S` is a no-op.

### Input

- `X` -- set
- `Type{S}` -- set type

### Output

The input set.
"""
function overapproximate(X::S, ::Type{S}) where {S <: LazySet}
    return X
end
