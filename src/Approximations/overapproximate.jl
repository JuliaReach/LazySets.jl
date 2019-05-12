import Base: convert

using LazySets: block_to_dimension_indices,
                substitute_blocks,
                get_constrained_lowdimset

"""
    overapproximate(X::S, ::Type{S}, args...) where {S<:LazySet}

Overapproximating a set of type `S` with type `S` is a no-op.

### Input

- `X`       -- set
- `Type{S}` -- target set type
- `args`    -- further arguments (ignored)

### Output

The input set.
"""
function overapproximate(X::S, ::Type{S}, args...) where {S<:LazySet}
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

- `S`        -- convex set, assumed to be two-dimensional
- `HPolygon` -- type for dispatch
- `ε`        -- (optional, default: `Inf`) error bound

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

Return a tight overapproximation of the Cartesian product array of a finite
number of convex sets with and hyperrectangle.

### Input

- `S`              -- Cartesian product array of a finite number of convex set
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

This method falls back to the corresponding `convert` method. Since the sets wrapped
by the Cartesian product array are hyperrectangles, it can be done efficiently
without overapproximation.
"""
function overapproximate(S::CartesianProductArray{N, <:AbstractHyperrectangle{N}},
                          ::Type{<:Hyperrectangle}) where {N<:Real}
    return convert(Hyperrectangle, S)
end

"""
    overapproximate(S::CartesianProduct{N, <:AbstractHyperrectangle{N}, <:AbstractHyperrectangle{N}},
                    ::Type{<:Hyperrectangle}) where {N<:Real}

Return a tight overapproximation of the Cartesian product of two
hyperrectangles by a new hyperrectangle.

### Input

- `S`              -- Cartesian product of two hyperrectangular sets
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

This method falls back to the corresponding `convert` method. Since the sets wrapped
by the Cartesian product are hyperrectangles, it can be done efficiently
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
    overapproximate(X::LazySet{N}, dir::AbstractDirections{N}) where {N}

Overapproximating a set with template directions.

### Input

- `X`   -- set
- `dir` -- (concrete) direction representation

### Output

A polyhedron overapproximating the set `X` with the directions from `dir`.
If the directions are known to be bounded, the result is an `HPolytope`,
otherwise the result is an `HPolyhedron`.
"""
function overapproximate(X::LazySet{N}, dir::AbstractDirections{N}) where {N}
    halfspaces = Vector{LinearConstraint{N}}()
    sizehint!(halfspaces, length(dir))
    T = isbounding(dir) ? HPolytope : HPolyhedron
    H = T(halfspaces)
    for d in dir
        addconstraint!(H, LinearConstraint(d, ρ(d, X)))
    end
    return H
end

"""
    overapproximate(X::LazySet{N}, dir::Type{<:AbstractDirections}) where {N}

Overapproximating a set with template directions.

### Input

- `X`   -- set
- `dir` -- type of direction representation

### Output

A polyhedron overapproximating the set `X` with the directions from `dir`.
If the directions are known to be bounded, the result is an `HPolytope`,
otherwise the result is an `HPolyhedron`.
"""
function overapproximate(X::LazySet{N},
                         dir::Type{<:AbstractDirections}) where {N}
    return overapproximate(X, dir{N}(dim(X)))
end

"""
    overapproximate(S::LazySet{N}, ::Type{<:Interval}) where {N<:Real}

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
function overapproximate(S::LazySet{N}, ::Type{<:Interval}) where {N<:Real}
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
                                           <:AbstractPolyhedron{N},
                                           <:AbstractPolyhedron{N}},
                         dir::AbstractDirections{N};
                         kwargs...
                        ) where {N<:Real}
    # important: the result may not be a polytope! 
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
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

# ==========================================
# Functionality that requires TaylorModels
# ==========================================

# function to be loaded by Requires
function load_taylormodels_overapproximation()

return quote

using .TaylorModels: Taylor1, TaylorN, TaylorModelN, TaylorModel1,
                     polynomial, remainder, domain,
                     normalize_taylor, linear_polynomial,
                     constant_term, evaluate, mid

# helper functions
@inline get_linear_coeffs(p::Taylor1) = linear_polynomial(p).coeffs[2:end]
@inline get_linear_coeffs(p::TaylorN) = linear_polynomial(p).coeffs[2].coeffs

"""
    overapproximate(vTM::Vector{TaylorModel1{T, S}},
                    ::Type{Zonotope}) where {T, S}

Overapproximate a taylor model in one variable with a zonotope.

### Input

- `vTM`      -- `TaylorModel1`
- `Zonotope` --  type for dispatch

### Output

A zonotope that overapproximates the range of the given taylor model.

### Examples

If the polynomials are linear, this functions exactly transforms to a zonotope.
However, the nonlinear case necessarily introduces overapproximation error.
Consider the linear case first:

```julia
julia> using LazySets, TaylorModels

julia> const IA = IntervalArithmetic;

julia> I = IA.Interval(-0.5, 0.5) # interval remainder
[-0.5, 0.5]

julia> x₀ = IA.Interval(0.0) # expansion point
[0, 0]

julia> D = IA.Interval(-3.0, 1.0)
[-3, 1]

julia> p1 = Taylor1([2.0, 1.0], 2) # define a linear polynomial
 2.0 + 1.0 t + 𝒪(t³)

julia> p2 = Taylor1([0.9, 3.0], 2) # define another linear polynomial
 0.9 + 3.0 t + 𝒪(t³)

julia> vTM = [TaylorModel1(pi, I, x₀, D) for pi in [p1, p2]]
2-element Array{TaylorModel1{Float64,Float64},1}:
 2.0 + 1.0 t + [-0.5, 0.5]
 0.9 + 3.0 t + [-0.5, 0.5]
```

Here, `vTM` is a taylor model vector, since each component is a taylor model in
one variable (`TaylorModel1`). Using `overapproximate(vTM, Zonotope)` we can
compute its associated zonotope in generator representation:

```julia
julia> using LazySets.Approximations

julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
2-element Array{Float64,1}:
  1.0
 -2.1

julia> Matrix(genmat(Z))
2×3 Array{Float64,2}:
 2.0  0.5  0.0
 6.0  0.0  0.5
```

Note how the generators of this zonotope mainly consist of two pieces: one comes
from the linear part of the polynomials, and another one that corresponds to the
interval remainder. This conversion gives the same upper and lower bounds as the
range evaluation using interval arithmetic:

```julia
julia> X = box_approximation(Z)
Hyperrectangle{Float64}([1.0, -2.1], [2.5, 6.5])

julia> Y = evaluate(vTM[1], vTM[1].dom) × evaluate(vTM[2], vTM[2].dom)
[-1.5, 3.5] × [-8.60001, 4.40001]

julia> H = convert(Hyperrectangle, Y) # this IntevalBox is the same as X
Hyperrectangle{Float64}([1.0, -2.1], [2.5, 6.5])
```
However, the zonotope returns better results if we want to approximate the `TM`,
since it is not axis-aligned:

```julia
julia> d = [-0.35, 0.93];

julia> ρ(d, Z) < ρ(d, X)
true
```

This function also works if the polynomials are non-linear; for example suppose
that we add a third polynomial with a quadratic term:

```julia
julia> p3 = Taylor1([0.9, 3.0, 1.0], 3);

julia> vTM = [TaylorModel1(pi, I, x₀, D) for pi in [p1, p2, p3]]
3-element Array{TaylorModel1{Float64,Float64},1}:
           2.0 + 1.0 t + [-0.5, 0.5]
           0.9 + 3.0 t + [-0.5, 0.5]
  0.9 + 3.0 t + 1.0 t² + [-0.5, 0.5]

julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
3-element Array{Float64,1}:
  1.0
 -2.1
  0.8999999999999999

julia> Matrix(genmat(Z))
3×4 Array{Float64,2}:
 2.0  0.5  0.0  0.0
 6.0  0.0  0.5  0.0
 6.0  0.0  0.0  6.5
```
The fourth and last generator corresponds to the addition between the interval
remainder and the box overapproximation of the nonlinear part of `p3` over the
domain.

### Algorithm

Let ``\\text{vTM} = (p, I)`` be a vector of ``m`` taylor models, where ``I``
is the interval remainder in ``\\mathbb{R}^m``. Let ``p_{lin}``
(resp. ``p_{nonlin}``) correspond to the linear (resp. nonlinear) part of each
scalar polynomial.

The range of ``\\text{vTM}`` can be enclosed by a zonotope with center ``c``
and matrix of generators ``G``, ``Z = ⟨c, G⟩``, by performing a conservative
linearization of ``\\text{vTM}``:

```math
    vTM' = (p', I') := (p_{lin} − p_{nonlin} , I + \\text{Int}(p_{nonlin})).
```

This algorithm proceeds in two steps:

1- Conservatively linearize ``\\text{vTM}`` as above and compute a box
   overapproximation of the nonlinear part.
2- Transform the linear taylor model to a zonotope exactly through variable
   normalization onto the symmetric intervals ``[-1, 1]``.
"""
function overapproximate(vTM::Vector{TaylorModel1{T, S}},
                            ::Type{Zonotope}) where {T, S}
    m = length(vTM)

    # preallocations
    c = Vector{T}(undef, m) # center of the zonotope
    gen_lin = Matrix{T}(undef, m, 1) # generator of the linear part
    gen_rem = Vector{T}(undef, m) # generators of the remainder

    # compute overapproximation
    return _overapproximate_vTM_zonotope!(vTM, c, gen_lin, gen_rem)
end

"""
    overapproximate(vTM::Vector{TaylorModelN{N, T, S}},
                    ::Type{Zonotope}) where {N,T, S}


Overapproximate a multivariate taylor model with a zonotope.

### Input

- `vTM`      -- `TaylorModelN`
- `Zonotope` -- type for dispatch

### Output

A zonotope that overapproximates the range of the given taylor model.

### Examples

Consider a vector of two 2-dimensional taylor models of order 2 and 4
respectively.

```julia
julia> using LazySets, LazySets.Approximations, TaylorModels

julia> const IA = IntervalArithmetic;

julia> x₁, x₂ = set_variables(Float64, ["x₁", "x₂"], order=8)
2-element Array{TaylorN{Float64},1}:
  1.0 x₁ + 𝒪(‖x‖⁹)
  1.0 x₂ + 𝒪(‖x‖⁹)

julia> x₀ = IntervalBox(0..0, 2) # expansion point
[0, 0] × [0, 0]

julia> Dx₁ = IA.Interval(0.0, 3.0) # domain for x₁
[0, 3]

julia> Dx₂ = IA.Interval(-1.0, 1.0) # domain for x₂
[-1, 1]

julia> D = Dx₁ × Dx₂ # take the Cartesian product of the domain on each variable
[0, 3] × [-1, 1]

julia> r = IA.Interval(-0.5, 0.5) # interval remainder
[-0.5, 0.5]

julia> p1 = 1 + x₁^2 - x₂
 1.0 - 1.0 x₂ + 1.0 x₁² + 𝒪(‖x‖⁹)

julia> p2 = x₂^3 + 3x₁^4 + x₁ + 1
 1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + 𝒪(‖x‖⁹)

julia> vTM = [TaylorModelN(pi, r, x₀, D) for pi in [p1, p2]]
2-element Array{TaylorModelN{2,Float64,Float64},1}:
             1.0 - 1.0 x₂ + 1.0 x₁² + [-0.5, 0.5]
   1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + [-0.5, 0.5]

julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
2-element Array{Float64,1}:
   5.5
 124.0

julia> Matrix(genmat(Z))
2×4 Array{Float64,2}:
 0.0  -1.0  5.0    0.0
 1.5   0.0  0.0  123.0
```

### Algorithm

We refer to the algorithm description for the univariate case.
"""
function overapproximate(vTM::Vector{TaylorModelN{N, T, S}},
                         ::Type{Zonotope}) where {N,T, S}
    m = length(vTM)
    n = N # number of variables is get_numvars() in TaylorSeries

    # preallocations
    c = Vector{T}(undef, m) # center of the zonotope
    gen_lin = Matrix{T}(undef, m, n) # generator of the linear part
    gen_rem = Vector{T}(undef, m) # generators for the remainder

    # compute overapproximation
    return _overapproximate_vTM_zonotope!(vTM, c, gen_lin, gen_rem)
end

function _overapproximate_vTM_zonotope!(vTM, c, gen_lin, gen_rem)
    @inbounds for (i, x) in enumerate(vTM)
        xpol, xdom = polynomial(x), domain(x)

        # linearize the TM
        pol_lin = constant_term(xpol) + linear_polynomial(xpol)
        rem_nonlin = remainder(x)

        # build an overapproximation of the nonlinear terms
        pol_nonlin = xpol - pol_lin
        rem_nonlin += evaluate(pol_nonlin, xdom)

        # normalize the linear polynomial to the symmetric interval [-1, 1]
        Q = normalize_taylor(pol_lin, xdom, true)

        # build the generators
        α = mid(rem_nonlin)
        c[i] = constant_term(Q) + α  # constant terms
        gen_lin[i, :] = get_linear_coeffs(Q) # linear terms
        gen_rem[i] = abs(rem_nonlin.hi - α)
    end
    return Zonotope(c, hcat(gen_lin, Diagonal(gen_rem)))
end

end # quote
end # load_taylormodels_overapproximation

# =============================================
# Functionality that requires IntervalMatrices
# =============================================

function load_intervalmatrices_overapproximation()
return quote

using .IntervalMatrices: AbstractIntervalMatrix, split

# temporary patch for IntervalArithmetic#317
function convert(::Type{IntervalMatrices.Interval{T}},
                 x::IntervalMatrices.Interval{T}) where {T<:Real}
    return x
end

"""
    overapproximate(lm::LinearMap{N, <:AbstractZonotope{N}, NM,
                                  <:AbstractIntervalMatrix{<:NM}},
                    ::Type{<:Zonotope})::Zonotope{N} where {N<:Real, NM}

Overapproximate an interval-matrix linear map of a zonotopic set by a new
zonotope.

### Input

- `lm`       -- interval-matrix linear map of a zonotopic set
- `Zonotope` -- type for dispatch

### Output

A zonotope overapproximating the linear map.

### Algorithm

This function implements the method proposed in [1].

Given an interval matrix ``M = \\tilde{M} + ⟨-\\hat{M},\\hat{M}⟩`` (split into a
conventional matrix and a symmetric interval matrix) and a zonotope
``⟨c, g_1, …, g_m⟩``, we compute the resulting zonotope
``⟨\\tilde{M}c, \\tilde{M}g_1, …, \\tilde{M}g_m, v_1, …, v_n⟩`` where the
``v_j``, ``j = 1, …, n``, are defined as

```math
    v_j = \\begin{cases} 0 & i ≠ j \\\\
          \\hat{M}_j (|c| + \\sum_{k=1}^m |g_k|) & i = j. \\end{cases}
```

[1] Althoff, Stursberg, Buss. Reachability analysis of linear systems with
uncertain parameters and inputs. CDC 2007.
"""
function overapproximate(lm::LinearMap{N, <:AbstractZonotope{N}, NM,
                                       <:AbstractIntervalMatrix{<:NM}},
                         ::Type{<:Zonotope})::Zonotope{N} where {N<:Real, NM}
    Mc, Ms = split(lm.M)
    Z = lm.X
    c = Mc * center(Z)
    n = dim(lm)
    nG = ngens(Z)
    G = zeros(N, n, nG + n)
    vector_sum = abs.(center(Z))
    @inbounds for (j, g) in enumerate(generators(Z))
        G[:, j] = Mc * g
        vector_sum += abs.(g)
    end
    @inbounds for i in 1:n
        row = @view Ms[i, :]
        G[i, i + nG] = dot(row, vector_sum)
    end
    return Zonotope(c, G)
end

end end  # quote / load_intervalmatrices_overapproximation()

# ==========================================
# Lazy linear maps of Cartesian products
# ==========================================

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                    ::Type{CartesianProductArray{N, S}}
                   ) where {N, S<:LazySet{N}}

Decompose a lazy linear map of a Cartesian product array while keeping the
original block structure.

### Input

- `lm`                    -- lazy linear map of Cartesian product array
- `CartesianProductArray` -- type for dispatch

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                         ::Type{CartesianProductArray{N, S}}
                        ) where {N, S<:LazySet{N}}
    cpa = array(lm.X)
    arr = Vector{S}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, S)
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                    ::Type{<:CartesianProductArray},
                    dir::Type{<:AbstractDirections}) where {N}

Decompose a lazy linear map of a Cartesian product array with template
directions while keeping the original block structure.

### Input

- `lm`                    -- lazy linear map of a Cartesian product array
- `CartesianProductArray` -- type for dispatch
- `dir`                   -- template directions for overapproximation

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                         ::Type{<:CartesianProductArray},
                         dir::Type{<:AbstractDirections}) where {N}
    cpa = array(lm.X)
    arr = Vector{HPolytope{N}}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, dir)
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                    ::Type{<:CartesianProductArray},
                    set_type::Type{<:LazySet}) where {N}

Decompose a lazy linear map of a Cartesian product array with a given set type
while keeping the original block structure.

### Input

- `lm`                    -- lazy linear map of a Cartesian product array
- `CartesianProductArray` -- type for dispatch
- `set_type`              -- set type for overapproximation

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                         ::Type{<:CartesianProductArray},
                         set_type::Type{<:LazySet}) where {N}
    cpa = array(lm.X)
    arr = Vector{set_type{N}}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, set_type)
end

function _overapproximate_lm_cpa!(arr, M, cpa, overapprox_option)
    # construct Minkowski sum for block row
    function _block_row(cpa::Vector{S}, M::AbstractMatrix{N},
                        row_range::UnitRange{Int}) where {N, S<:LazySet{N}}
        arr_inner = Vector{LinearMap{N, <:S}}(undef, length(cpa))
        col_start_ind, col_end_ind = 1, 0
        @inbounds for (j, bj) in enumerate(cpa)
            col_end_ind += dim(bj)
            arr_inner[j] =
                LinearMap(M[row_range, col_start_ind:col_end_ind], bj)
            col_start_ind = col_end_ind + 1
        end
        return MinkowskiSumArray(arr_inner)
    end

    row_start_ind, row_end_ind = 1, 0
    @inbounds for (i, bi) in enumerate(cpa)
        row_end_ind += dim(bi)
        ms = _block_row(cpa, M, row_start_ind:row_end_ind)
        arr[i] = overapproximate(ms, overapprox_option)
        row_start_ind = row_end_ind + 1
    end

    return CartesianProductArray(arr)
end

"""
    overapproximate(cap::Intersection{N,
                                      <:CartesianProductArray{N},
                                      <:AbstractPolyhedron{N}},
                    ::Type{CartesianProductArray}, oa) where {N}

Return the intersection of the Cartesian product of a finite number of convex
sets and a polyhedron.

### Input

 - `cap` -- Lazy intersection of Cartesian product array and polyhedron
 - `CartesianProductArray` -- type for dispatch
 - `oa`  -- overapproximation option

### Output

A `CartesianProductArray` that overapproximates the intersection of `cpa` and
`P`.

### Algorithm

The intersection only needs to be computed in the blocks of `cpa` that are
constrained in `P`.
Hence we first collect those constrained blocks in a lower-dimensional Cartesian
product array and then convert to an `HPolytope` `X`.
Then we take the intersection of `X` and the projection of `Y` onto the
corresponding dimensions.
(This projection is purely syntactic and exact.)
Finally we decompose the result again and plug together the unaffected old
blocks and the newly computed blocks.
The result is a `CartesianProductArray` with the same block structure as in `X`.
"""
function overapproximate(cap::Intersection{N,
                                           <:CartesianProductArray{N},
                                           <:AbstractPolyhedron{N}},
                            ::Type{CartesianProductArray}, oa) where {N}

    cpa, P = cap.X, cap.Y

    cpa_low_dim, vars, block_structure, blocks =
        get_constrained_lowdimset(cpa, P)

    hpoly_low_dim = HPolytope(constraints_list(cpa_low_dim))
    low_intersection = intersection(hpoly_low_dim, project(P, vars))

    if isempty(low_intersection)
        return EmptySet{N}()
    end

    decomposed_low_set = decompose(low_intersection, block_structure, oa)

    return substitute_blocks(decomposed_low_set, cpa, blocks)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                            <:AbstractPolyhedron{N},
                                            <:CartesianProductArray{N}},
                            ::Type{CartesianProductArray}, oa) where {N}
    overapproximate(Intersection(cap.Y, cap.X), oa)
end
