"""
    overapproximate(X::S, ::Type{S}) where {S<:LazySet}

Overapproximating a set of type `S` with type `S` is a no-op.

### Input

- `X` -- set
- `Type{S}` -- set type

### Output

The input set.
"""
function overapproximate(X::S, ::Type{S}) where {S<:LazySet}
    return X
end

"""
    overapproximate(S::LazySet{N},
                    ::Type{<:HPolygon},
                    [ε]::Real=Inf)::HPolygon where {N<:Real}

Return an approximation of a given 2D convex set.
If no error tolerance is given, or is `Inf`, the result is a box-shaped polygon.
Otherwise the result is an ε-close approximation as a polygon.

### Input

- `S`           -- convex set, assumed to be two-dimensional
- `HPolygon`    -- type for dispatch
- `ε`           -- (optional, default: `Inf`) error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(S::LazySet{N},
                         ::Type{<:HPolygon},
                         ε::Real=Inf
                        )::HPolygon where {N<:Real}
    @assert dim(S) == 2
    if ε == Inf
        pe = σ(DIR_EAST(N), S)
        pn = σ(DIR_NORTH(N), S)
        pw = σ(DIR_WEST(N), S)
        ps = σ(DIR_SOUTH(N), S)
        constraints = Vector{LinearConstraint{N}}(undef, 4)
        constraints[1] = LinearConstraint(DIR_EAST(N), dot(pe, DIR_EAST(N)))
        constraints[2] = LinearConstraint(DIR_NORTH(N), dot(pn, DIR_NORTH(N)))
        constraints[3] = LinearConstraint(DIR_WEST(N), dot(pw, DIR_WEST(N)))
        constraints[4] = LinearConstraint(DIR_SOUTH(N), dot(ps, DIR_SOUTH(N)))
        return HPolygon(constraints, sort_constraints=false)
    else
        return tohrep(approximate(S, ε))
    end
end

"""
    overapproximate(S::LazySet, ε::Real)::HPolygon

Alias for `overapproximate(S, HPolygon, ε)`.
"""
function overapproximate(S::LazySet,
                         ε::Real;
                        )::HPolygon
    return overapproximate(S, HPolygon, ε)
end

"""
    overapproximate(S::LazySet,
                    Type{<:Hyperrectangle})::Union{Hyperrectangle, EmptySet}

Return an approximation of a given set as a hyperrectangle.

### Input

- `S`              -- set
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.
"""
function overapproximate(S::LazySet,
                         ::Type{<:Hyperrectangle};
                        )::Union{Hyperrectangle, EmptySet}
    box_approximation(S)
end

"""
    overapproximate(S::LazySet)::Union{Hyperrectangle, EmptySet}

Alias for `overapproximate(S, Hyperrectangle)`.
"""
overapproximate(S::LazySet;
               )::Union{Hyperrectangle, EmptySet} =
    overapproximate(S, Hyperrectangle)

"""
    overapproximate(S::ConvexHull{N, Zonotope{N}, Zonotope{N}},
                    ::Type{<:Zonotope})::Zonotope where {N<:Real}

Overapproximate the convex hull of two zonotopes.

### Input

- `S`           -- convex hull of two zonotopes of the same order
- `Zonotope`    -- type for dispatch

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

"""
    overapproximate(Z::Zonotope, ::Type{<:Hyperrectangle})::Hyperrectangle

Return a tight overapproximation of a zonotope with an axis-aligned box.

### Input

- `Z`              -- zonotope
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

This function implements the method in [Section 5.1.2, 1]. A zonotope
``Z = ⟨c, G⟩`` can be overapproximated tightly by an axis-aligned box
(i.e. a `Hyperrectangle`) such that its center is ``c`` and the radius along
dimension ``i`` is the column-sum of the absolute values of the ``i``-th row
of ``G`` for ``i = 1,…, p``, where ``p`` is the number of generators of ``Z``.

[1] *Althoff, M., Stursberg, O., & Buss, M. (2010). Computing reachable sets of
hybrid systems using a combination of zonotopes and polytopes. Nonlinear analysis:
hybrid systems, 4(2), 233-249.*
"""
function overapproximate(Z::Zonotope, ::Type{<:Hyperrectangle})::Hyperrectangle
    r = sum(abs.(Z.generators), dims=2)[:]
    return Hyperrectangle(Z.center, r)
end

"""
    overapproximate(X::LazySet{N}, dir::AbstractDirections{N})::HPolytope{N}
        where {N}

Overapproximating a set with template directions.

### Input

- `X`           -- set
- `dir`         -- direction representation

### Output

A `HPolytope` overapproximating the set `X` with the directions from `dir`.
"""
function overapproximate(X::LazySet{N},
                         dir::AbstractDirections{N}
                        )::HPolytope{N} where {N}
    halfspaces = Vector{LinearConstraint{N}}()
    sizehint!(halfspaces, length(dir))
    H = HPolytope(halfspaces)
    for d in dir
        addconstraint!(H, LinearConstraint(d, ρ(d, X)))
    end
    return H
end

"""
    overapproximate(S::LazySet{N}, ::Type{Interval}) where {N<:Real}

Return the overapproximation of a real unidimensional set with an interval.

### Input

- `S`           -- one-dimensional set
- `Interval`    -- type for dispatch

### Output

An interval.
"""
function overapproximate(S::LazySet{N}, ::Type{Interval}) where {N<:Real}
    @assert dim(S) == 1
    lo = σ([-one(N)], S)[1]
    hi = σ([one(N)], S)[1]
    return Interval(lo, hi)
end

function overapproximate_cap_helper(X::LazySet{N},                # compact set
                                    P::Union{AbstractPolytope{N}, # polyhedron
                                             HPolyhedron{N}},
                                    dir::AbstractDirections{N};
                                    kwargs...
                                   ) where {N<:Real}
    Hi = constraints_list(P)
    m = length(Hi)
    Q = HPolytope{N}()

    for di in dir
        ρ_X_Hi_min = ρ(di, X ∩ Hi[1], kwargs...)
        for i in 2:m
            ρ_X_Hi = ρ(di, X ∩ Hi[i], kwargs...)
            if ρ_X_Hi < ρ_X_Hi_min
                ρ_X_Hi_min = ρ_X_Hi
            end
        end
        addconstraint!(Q, HalfSpace(di, ρ_X_Hi_min))
    end
    return Q
end

"""
    overapproximate(cap::Intersection{N,
                                      <:LazySet,
                                      <:Union{AbstractPolytope{N}, HPolyhedron{N}}},
                    dir::AbstractDirections{N};
                    kwargs...
                   ) where {N<:Real}

Return the overapproximation of the intersection between a compact set and a
polytope given a set of template directions.

### Input

- `cap`         -- intersection of a compact set and a polytope
- `dir`         -- template directions
- `kwargs`      -- additional arguments that are passed to the support function
                   algorithm

### Output

A polytope in H-representation such that the normal direction of each half-space
is given by an element of `dir`.

### Algorithm

Let `di` be a direction drawn from the set of template directions `dir`.
Let `X` be the compact set and let `P` be the polytope. We overapproximate the
set `X ∩ H` with a polytope in constraint representation using a given set of
template directions `dir`.

The idea is to solve the univariate optimization problem `ρ(di, X ∩ Hi)` for
each half-space in the set `P` and then take the minimum.
This gives an overapproximation of the exact support function.

This algorithm is inspired from [G. Frehse, R. Ray. Flowpipe-Guard Intersection
for Reachability Computations with Support
Functions](https://www.sciencedirect.com/science/article/pii/S1474667015371809).

### Notes

This method relies on having available the `constraints_list` of the polytope
`P`.

This method of overapproximations can return a non-empty set even if the original
intersection is empty.
"""
function overapproximate(cap::Intersection{N,
                                           <:LazySet,
                                           <:Union{AbstractPolytope{N}, HPolyhedron{N}}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:Union{AbstractPolytope{N}, HPolyhedron{N}},
                                           <:LazySet},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope{N},
                                           <:Union{AbstractPolytope{N}, HPolyhedron{N}}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:Union{AbstractPolytope{N}, HPolyhedron{N}},
                                           <:AbstractPolytope{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope{N},
                                           <:HPolyhedron{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:HPolyhedron{N},
                                           <:AbstractPolytope{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

"""
    overapproximate(cap::Intersection{N, <:HalfSpace{N}, <:AbstractPolytope{N}},
                    dir::AbstractDirections{N};
                    [kwargs]...
                   ) where {N<:Real}

Return the overapproximation of the intersection between a half-space and a
polytope given a set of template directions.

### Input

- `cap`         -- intersection of a half-space and a polytope
- `dir`         -- template directions
- `kwargs`      -- additional arguments that are passed to the support function
                   algorithm

### Output

A polytope in H-representation such that the normal direction of each half-space
is given by an element of `dir`.
"""
function overapproximate(cap::Intersection{N,
                                           <:HalfSpace{N},
                                           <:AbstractPolytope{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    H = HPolytope{N}()
    c = constraints_list(H)
    append!(c, constraints_list(cap.X))
    append!(c, constraints_list(cap.Y))
    return overapproximate(H, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope{N},
                                           <:HalfSpace{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate(swap(cap), dir; kwargs...)
end
