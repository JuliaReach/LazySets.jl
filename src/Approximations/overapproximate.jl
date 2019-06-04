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
        constraints = Vector{LinearConstraint{N}}(undef, 4)
        constraints[1] = LinearConstraint(DIR_EAST(N), ρ(DIR_EAST(N), S))
        constraints[2] = LinearConstraint(DIR_NORTH(N), ρ(DIR_NORTH(N), S))
        constraints[3] = LinearConstraint(DIR_WEST(N), ρ(DIR_WEST(N), S))
        constraints[4] = LinearConstraint(DIR_SOUTH(N), ρ(DIR_SOUTH(N), S))
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
    return box_approximation(S)
end

"""
    overapproximate(S::CartesianProductArray{N, <:AbstractHyperrectangle{N}},
                    ::Type{<:Hyperrectangle}) where {N<:Real}

Return a tight overapproximation of the cartesian product array of a finite
number of convex sets with and hyperrectangle.

### Input

- `S`              -- cartesian product array of a finite number of convex set
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

This method falls back to the corresponding `convert` method. Since the sets wrapped
by the cartesian product array are hyperrectangles, it can be done efficiently
without overapproximation.
"""
function overapproximate(S::CartesianProductArray{N, <:AbstractHyperrectangle{N}},
                          ::Type{<:Hyperrectangle}) where {N<:Real}
    return convert(Hyperrectangle, S)
end

"""
    overapproximate(S::CartesianProduct{N, <:AbstractHyperrectangle{N}, <:AbstractHyperrectangle{N}},
                    ::Type{<:Hyperrectangle}) where {N<:Real}

Return a tight overapproximation of the cartesian product of two
hyperrectangles by a new hyperrectangle.

### Input

- `S`              -- cartesian product of two hyperrectangular sets
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

This method falls back to the corresponding `convert` method. Since the sets wrapped
by the cartesian product are hyperrectangles, it can be done efficiently
without overapproximation.
"""
function overapproximate(S::CartesianProduct{N, <:AbstractHyperrectangle{N}, <:AbstractHyperrectangle{N}},
                          ::Type{<:Hyperrectangle}) where {N<:Real}
    return convert(Hyperrectangle, S)
end

"""
    overapproximate(lm::LinearMap{N, <:AbstractHyperrectangle{N}},
                    ::Type{Hyperrectangle}) where {N}

Return a tight overapproximation of the linear map of a hyperrectangular set
using a hyperrectangle.

### Input

- `S`              -- linear map of a hyperrectangular set
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

If `c` and `r` denote the center and vector radius of a hyperrectangle `H`,
a tight hyperrectangular overapproximation of `M * H` is obtained by transforming
`c ↦ M*c` and `r ↦ abs.(M) * c`, where `abs.(⋅)` denotes the element-wise absolute
value operator. 
"""
function overapproximate(lm::LinearMap{N, <:AbstractHyperrectangle{N}},
                         ::Type{Hyperrectangle}) where {N<:Real}
    M, X = lm.M, lm.X
    center_MX = M * center(X)
    radius_MX = abs.(M) * radius_hyperrectangle(X)
    return Hyperrectangle(center_MX, radius_MX)
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

- `S`        -- convex hull of two zonotopes
- `Zonotope` -- type for dispatch

### Algorithm

This function implements the method proposed in [1].
The convex hull of two zonotopes ``Z₁`` and ``Z₂`` of the same order,
that we write

```math
Z_j = ⟨c^{(j)}, g^{(j)}_1, …, g^{(j)}_p⟩
```
for ``j = 1, 2``, can be overapproximated as follows:

```math
CH(Z_1, Z_2) ⊆ \\frac{1}{2}⟨c^{(1)}+c^{(2)}, g^{(1)}_1+g^{(2)}_1, …, g^{(1)}_p+g^{(2)}_p, c^{(1)}-c^{(2)}, g^{(1)}_1-g^{(2)}_1, …, g^{(1)}_p-g^{(2)}_p⟩.
```

If the zonotope order is not the same, this algorithm calls
`reduce_order` to reduce the order to the minimum of the arguments.

It should be noted that the output zonotope is not necessarily the minimal
enclosing zonotope, which is in general expensive in high dimensions. This is
further investigated in [2].

[1] Reachability of Uncertain Linear Systems Using Zonotopes, A. Girard.
    HSCC 2005.

[2] Zonotopes as bounding volumes, L. J. Guibas et al, Proc. of Symposium on
    Discrete Algorithms, pp. 803-812.
"""
function overapproximate(S::ConvexHull{N, Zonotope{N}, Zonotope{N}},
                         ::Type{<:Zonotope})::Zonotope where {N<:Real}
    Z1, Z2 = S.X, S.Y
    
    # reduce to the same order if possible
    if order(Z1) != order(Z2)
        min_order = min(order(Z1), order(Z2))
        Z1 = reduce_order(Z1, min_order)
        Z2 = reduce_order(Z2, min_order)
    end

    if order(Z1) >= order(Z2)
        c, G = _overapproximate_convex_hull_zonotope(Z1, Z2)
    else
        c, G = _overapproximate_convex_hull_zonotope(Z2, Z1)
    end
    return Zonotope(c, G)
end

# assumes that dim(Z1) == dim(Z2) and order(Z1) >= order(Z2)
function _overapproximate_convex_hull_zonotope(Z1::Zonotope{N}, Z2::Zonotope{N}) where {N}
    c = (Z1.center + Z2.center)/N(2)

    # the case of equal order is treated separately to avoid a slicing (this creates a copy)
    if order(Z1) == order(Z2)        
        G = hcat(Z1.generators .+ Z2.generators,
                 Z1.center - Z2.center,
                 Z1.generators .- Z2.generators)/N(2)
    else
        G = hcat((Z1.generators[:, 1:ngens(Z2)] .+ Z2.generators)/N(2),
                 (Z1.center - Z2.center)/N(2),
                 (Z1.generators[:, 1:ngens(Z2)] .- Z2.generators)/N(2),
                 Z1.generators[:, ngens(Z2)+1:end])
    end
    return c, G
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

- `X`   -- set
- `dir` -- (concrete) direction representation

### Output

An `HPolytope` overapproximating the set `X` with the directions from `dir`.
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
    overapproximate(X::LazySet{N},
                    dir::Type{<:AbstractDirections})::HPolytope{N} where {N}

Overapproximating a set with template directions.

### Input

- `X`   -- set
- `dir` -- type of direction representation

### Output

A `HPolytope` overapproximating the set `X` with the directions from `dir`.
"""
function overapproximate(X::LazySet{N},
                         dir::Type{<:AbstractDirections}
                        )::HPolytope{N} where {N}
    return overapproximate(X, dir{N}(dim(X)))
end

"""
    overapproximate(S::LazySet{N}, ::Type{Interval}) where {N<:Real}

Return the overapproximation of a real unidimensional set with an interval.

### Input

- `S`        -- one-dimensional set
- `Interval` -- type for dispatch

### Output

An interval.

### Algorithm

The method relies on the exact conversion to `Interval`. Two support
function evaluations are needed in general.
"""
function overapproximate(S::LazySet{N}, ::Type{Interval}) where {N<:Real}
    @assert dim(S) == 1 "cannot overapproximate a $(dim(S))-dimensional set with an `Interval`"
    return convert(Interval, S)
end

function overapproximate_cap_helper(X::LazySet{N},             # compact set
                                    P::AbstractPolyhedron{N},  # polyhedron
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
    overapproximate(cap::Intersection{N, <:LazySet, <:AbstractPolyhedron{N}},
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
                                           <:AbstractPolyhedron{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron{N},
                                           <:LazySet},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope{N},
                                           <:AbstractPolyhedron{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron{N},
                                           <:AbstractPolytope{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope{N},
                                           <:AbstractPolytope{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
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

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray{N,
                            <:LazySet{N}}}, ::Type{<:CartesianProductArray}, oa=Hyperrectangle)  where {N}

Overapproximating a lazy linear map of cartesian product array with template directions for each block.

### Input

- `lm`   -- lazy linear map of cartesian product array
- Type{<:CartesianProductArray} -- type to specify decomposed overapproximation
- `oa` -- (concrete) direction representation

### Output

An `CartesianProductArray` with overapproximation of the each block with the directions from `dir`.
"""
function overapproximate(lm::LinearMap{N, <:CartesianProductArray{N,
                            <:LazySet{N}}}, ::Type{<:CartesianProductArray}, oa=Hyperrectangle) where {N}

    @assert size(lm.M, 2) == dim(lm.X) "matrix needs to be commensurate with the cartesian product"

    cp = lm.X
    M = lm.M

    col_end_ind = 0

    return decomposed_overapproximation(cp, M, col_end_ind, oa)
end

#template directions
function decomposed_overapproximation(cpa::CartesianProductArray{N,<:LazySet{N}},
                                      M::AbstractMatrix{N},
                                      col_end_ind::Int, oa) where {N}
    array = allocate_result(oa, N)
    sizehint!(array,length(cpa.array))

    for i in 1:size(M, 2)
        ms = overapproximate_row_blocks(cpa, M, i)
        push!(array, overapproximate(ms, oa))
    end

    result = CartesianProductArray(array)
    return result
end

#overapproximation of linear map to each block
function overapproximate_row_blocks(cpa::CartesianProductArray{N,<:LazySet{N}},
                                    M::AbstractMatrix{N}, i::Int) where {N}
    col_start_ind, col_end_ind = 1, 0
    h_min_sum = MinkowskiSumArray()
    for bi in cpa.array
        col_end_ind += dim(bi)
        push!(h_min_sum.array, LinearMap(M[i : i, col_start_ind : col_end_ind], bi))
        col_start_ind = col_end_ind + 1
    end

    return h_min_sum
end

function allocate_result(oa, N)
    return Vector{LazySet{N}}()
end

function allocate_result(oa::Type{<:LazySet}, N)
    return Vector{oa{N}}()
end

function allocate_result(oa::Type{<:AbstractDirections}, N)
    return Vector{HPolytope{N}}()
end
