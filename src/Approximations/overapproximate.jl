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



function load_taylormodels_overapproximation()  # function to be loaded by Requires
return quote
using .TaylorModels:Taylor1,TaylorN,TaylorModelN,TaylorModel1,normalize_taylor
using .TaylorModels:linear_polynomial,constant_term,evaluate,mid

"""
         overapproximate(vTM::Vector{TaylorModel1{T, S}},
                            ::Type{Zonotope}) where {T, S}

Overapproximate a taylor model in one variable with a zonotope.

### Input


- `vTM` -- `TaylorModel1`
- `Zonotope` -- type for dispatch

### Output

A zonotope that overapproximates the range of the given taylor model.

### Examples

```julia
julia> using TaylorModels, IntervalArithmetic


julia> δ = 0.5; I = Interval(-δ, δ)
[-0.5, 0.5]

julia> x₀ = Interval(0.0)
[0, 0]

julia> D = Interval(-3.0, 1.0)
[-3, 1]

julia> p = Taylor1([2.0, 1.0], 2)
 2.0 + 1.0 t + 𝒪(t³)

julia> p1 = Taylor1([0.9, 3.0], 2)
 0.9 + 3.0 t + 𝒪(t³)

julia> TM = [TaylorModel1(p, I, x₀, D),TaylorModel1(p1, I, x₀, D)]
2-element Array{TaylorModel1{Float64,Float64},1}:
 2.0 + 1.0 t + [-0.5, 0.5]
 0.9 + 3.0 t + [-0.5, 0.5]

julia> overapproximate(TM,Zonotope)
Zonotope{Float64}([1.0, -2.1],
   [1, 1]  =  2.0
   [2, 1]  =  6.0
   [1, 2]  =  0.5
   [2, 3]  =  0.5)
```

### Algorithm

The TaylorModel ``TM`` can be enclosed by a zonotope ``Z = ⟨c, G⟩``.
A simple but effective way to do so is to perform a conservative linearization
on the ``TM``,
``TM = (p′, I′) = (p − pN , I + Int(pN ))``
Now idea is to normalize the linear part and take Box overapproximate of nonlinear
part then collect center and generators of zonotope overapproximating TaylorModel.

"""
function overapproximate(vTM::Vector{TaylorModel1{T, S}},
                            ::Type{Zonotope}) where {T, S}
    m = length(vTM)
    c = Vector{T}(undef, m) # center of the zonotope
    gen = Vector{T}(undef, m) # generator of the linear part
    rem_gen = Vector{T}(undef, m) # generators for the remainder

    @inbounds for (i, x) in enumerate(vTM)
        # linearize the TM
        pol_lin = constant_term(x.pol) + linear_polynomial(x.pol)
        rem_nonlin = x.rem
        # if there are nonlinear terms, build an overapproximation
        if length(x.pol.coeffs) > 2
            pol_nonlin = x.pol - pol_lin
            rem_nonlin += evaluate(pol_nonlin, x.dom)
        end
        # normalize the linear polynomial to the symmetric interval [-1, 1]
        Q = normalize_taylor(pol_lin, x.dom, true)
        # build the generators
        α = mid(rem_nonlin)
        c[i] = Q.coeffs[1] + α
        gen[i] = Q.coeffs[2]
        rem_gen[i] = abs(rem_nonlin.hi - α)
    end
    return Zonotope(c, hcat(gen, Diagonal(rem_gen)))
end

"""
          overapproximate(vTM::Vector{TaylorModelN{N, T, S}},
                             ::Type{Zonotope}) where {N,T, S}


Overapproximate a multivariable taylor model with a zonotope.

### Input


- `vTM` -- `TaylorModelN`
- `Zonotope` -- type for dispatch

### Output

A zonotope that overapproximates the range of the given taylor model.

### Examples


```julia
julia> m = 4
4

julia> x₁, x₂ = set_variables(Float64, ["x₁", "x₂"], order=2*m)
2-element Array{TaylorN{Float64},1}:
  1.0 x₁ + 𝒪(‖x‖⁹)
  1.0 x₂ + 𝒪(‖x‖⁹)

julia> x₀ = Interval(0.0, 0.0) × Interval(0.0, 0.0)
[0, 0] × [0, 0]

julia> Dx₁ = Interval(0.0, 3.0)
[0, 3]

julia> Dx₂ = Interval(-1.0, 1.0)
[-1, 1]

julia> D = Dx₁ × Dx₂
[0, 3] × [-1, 1]

julia> δ = 0.5; I = Interval(-δ, δ)
[-0.5, 0.5]

julia> p = 1 + x₁^2 - x₂
 1.0 - 1.0 x₂ + 1.0 x₁² + 𝒪(‖x‖⁹)

julia> p1 = x₂^3 + 3x₁^4 + x₁ + 1
 1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + 𝒪(‖x‖⁹)

julia> TM = [ TaylorModelN(p, I, x₀, D), TaylorModelN(p1, I, x₀, D) ]
2-element Array{TaylorModelN{2,Float64,Float64},1}:
             1.0 - 1.0 x₂ + 1.0 x₁² + [-0.5, 0.5]
   1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + [-0.5, 0.5]

julia> overapproximate(TM, Zonotope)
Zonotope{Float64}([5.5, 124.0],
   [2, 1]  =  1.5
   [1, 2]  =  -1.0
   [1, 3]  =  5.0
   [2, 4]  =  123.0)
 ```
### Algorithm

The TaylorModel ``TM`` can be enclosed by a zonotope ``Z = ⟨c, G⟩``.
A simple but effective way to do so is to perform a conservative linearization
on the ``TM``,
``TM = (p′, I′) = (p − pN , I + Int(pN ))``
Now idea is to normalize the linear part and take Box overapproximate of nonlinear
part then collect center and generators of zonotope overapproximating TaylorModel.

"""
function overapproximate(vTM::Vector{TaylorModelN{N, T, S}},
                            ::Type{Zonotope}) where {N,T, S}
    m = length(vTM)
    n = length(vTM[1].pol.coeffs[2].coeffs)
    c = Vector{T}(undef, m) # center of the zonotope
    gen = Matrix{T}(undef, m, n) # generator of the linear part
    rem_gen = Vector{T}(undef, m) # generators for the remainder

    @inbounds for (i, x) in enumerate(vTM)
        # linearize the TM
        pol_lin = constant_term(x.pol) + linear_polynomial(x.pol)
        rem_nonlin = x.rem

        pol_nonlin = x.pol - pol_lin
        rem_nonlin += evaluate(pol_nonlin, x.dom)
        # normalize the linear polynomial to the symmetric interval [-1, 1]
        Q = normalize_taylor(pol_lin, x.dom, true)
        # build the generators
        α = mid(rem_nonlin)
        c[i] = Q.coeffs[1].coeffs[1] + α
        gen[i,1:n] = Q.coeffs[2].coeffs
        rem_gen[i] = abs(rem_nonlin.hi - α)
    end
    return Zonotope(c, hcat(gen, Diagonal(rem_gen)))
end
end
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray{N}},
                    ::Type{CartesianProductArray{N, S}}
                   ) where {N, S<:LazySet{N}}

Decompose a lazy linear map of a cartesian product array while keeping the
original block structure.

### Input

- `lm` -- lazy linear map of cartesian product array
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

Decompose a lazy linear map of a cartesian product array with template
directions while keeping the original block structure.

### Input

- `lm`  -- lazy linear map of a cartesian product array
- `CartesianProductArray` -- type for dispatch
- `dir` -- template directions for overapproximation

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

Decompose a lazy linear map of a cartesian product array with a given set type
while keeping the original block structure.

### Input

- `lm`  -- lazy linear map of a cartesian product array
- `CartesianProductArray` -- type for dispatch
- `set_type` -- set type for overapproximation

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
