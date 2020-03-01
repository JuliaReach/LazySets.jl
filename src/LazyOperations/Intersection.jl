import Base: isempty, ∈, ∩

export Intersection,
       isempty_known,
       set_isempty!,
       swap,
       use_precise_ρ,
       constraints_list,
       linear_map

"""
    IntersectionCache

Container for information cached by a lazy `Intersection` object.

### Fields

- `isempty` -- is the intersection empty? There are three possible states,
               encoded as `Int8` values -1, 0, 1:

    * ``-1`` - it is currently unknown whether the intersection is empty or not
    *  ``0`` - intersection is not empty
    *  ``1`` - intersection is empty
"""
mutable struct IntersectionCache
    isempty::Int8

    # default constructor
    IntersectionCache() = new(Int8(-1))
end

function isempty_known(c::IntersectionCache)
    return c.isempty != Int8(-1)
end

function isempty(c::IntersectionCache)
    @assert isempty_known(c) "'isempty_known' only works if 'isempty' " *
        "returns 'true'"
    return c.isempty == Int8(1)
end

function set_isempty!(c::IntersectionCache, isempty::Bool)
    c.isempty = isempty ? Int8(1) : Int8(0)
end

"""
    Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two convex sets.

### Fields

- `X`     -- convex set
- `Y`     -- convex set
- `cache` -- internal cache for avoiding recomputation; see
             [`IntersectionCache`](@ref)

### Examples

Create an expression, ``Z``, which lazily represents the intersection of two
squares ``X`` and ``Y``:

```jldoctest lazy_intersection
julia> X, Y = BallInf([0,0.], 0.5), BallInf([1,0.], 0.75);

julia> Z = X ∩ Y;

julia> typeof(Z)
Intersection{Float64,BallInf{Float64,Array{Float64,1}},BallInf{Float64,Array{Float64,1}}}

julia> dim(Z)
2
```

We can check if the intersection is empty with `isempty`:

```jldoctest lazy_intersection
julia> isempty(Z)
false
```

Do not confuse `Intersection` with the concrete operation, which is computed
with the lowercase `intersection` function:

```jldoctest lazy_intersection
julia> W = intersection(X, Y)
Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}([0.375, 0.0], [0.125, 0.5])
```
"""
struct Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
    cache::IntersectionCache

    # default constructor with dimension check
    function Intersection(X::S1, Y::S2;
                          cache::IntersectionCache=IntersectionCache()
                         ) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in an intersection must have the same " *
            "dimension"
        return new{N, S1, S2}(X, Y, cache)
    end
end

isoperationtype(::Type{<:Intersection}) = true
isconvextype(::Type{Intersection{N, S1, S2}}) where {N, S1, S2} = isconvextype(S1) && isconvextype(S2)

# Universe is the neutral element for Intersection
@neutral(Intersection, Universe)

# EmptySet is the absorbing element for Intersection
@absorbing(Intersection, EmptySet)

"""
    ∩

Alias for `Intersection`.
"""
∩(X::LazySet, Y::LazySet) = Intersection(X, Y)


# --- cache propagation functions ---


"""
    isempty_known(cap::Intersection)

Ask whether the status of emptiness is known.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` iff the emptiness status is known.
In this case, `isempty(cap)` can be used to obtain the status.
"""
function isempty_known(cap::Intersection)
    return isempty_known(cap.cache)
end

"""
    set_isempty!(cap::Intersection, isempty::Bool)

Set the status of emptiness in the cache.

### Input

- `cap`     -- intersection of two convex sets
- `isempty` -- new status of emptiness
"""
function set_isempty!(cap::Intersection, isempty::Bool)
    return set_isempty!(cap.cache, isempty)
end

"""
    swap(cap::Intersection{N, S1, S2}) where {N<:Real, S1, S2}

Return a new `Intersection` object with the arguments swapped.

### Input

- `cap` -- intersection of two convex sets

### Output

A new `Intersection` object with the arguments swapped.
The old cache is shared between the old and new objects.

### Notes

The advantage of using this function instead of manually swapping the arguments
is that the cache is shared.
"""
function swap(cap::Intersection{N, S1, S2}) where {N<:Real, S1, S2}
    return Intersection(cap.Y, cap.X, cache=cap.cache)
end


# --- LazySet interface functions ---


"""
    dim(cap::Intersection)

Return the dimension of an intersection of two convex sets.

### Input

- `cap` -- intersection of two convex sets

### Output

The ambient dimension of the intersection of two convex sets.
"""
function dim(cap::Intersection)
    return dim(cap.X)
end

"""
    σ(d::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}

Return the support vector of an intersection of two convex sets in a given
direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two convex sets

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}
    error("the exact support vector of an intersection is not implemented")
end

"""
    ρ(d::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}

Return an upper bound on the support function of the intersection of two convex
sets in a given direction.

### Input

- `d`    -- direction
- `cap`  -- intersection of two convex sets

### Output

An uper bound on the support function in the given direction.

### Algorithm

The support function of an intersection of ``X`` and ``Y`` is upper bounded by
the minimum of the support functions of ``X`` and ``Y``.
"""
function ρ(d::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}
    return min(ρ(d, cap.X), ρ(d, cap.Y))
end

function ρ_helper(d::AbstractVector{N},
                  cap::Intersection{N,
                                    <:LazySet{N},
                                    <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}},
                  algorithm::String;
                  kwargs...) where {N<:Real}
    @assert isbounded(cap.X) "the first set in the intersection must be bounded"
    X = cap.X # compact set
    H = cap.Y # halfspace or hyperplane or line

    # if the intersection is empty => stop
    if isempty(cap)
        error("the intersection is empty")
    end

    if !use_precise_ρ(cap) || algorithm == "simple"
        return invoke(ρ, Tuple{typeof(d), Intersection{N}}, d, cap)
    elseif algorithm == "line_search"
        require(:Optim; fun_name="ρ", explanation="(algorithm $algorithm)")
        (s, _) = _line_search(d, X, H; kwargs...)
        return s
    elseif algorithm == "projection"
        @assert H isa Hyperplane "the algorithm $algorithm cannot be used " *
            "with a $(typeof(H)); it only works with hyperplanes"
        return _projection(d, X, H; kwargs...)
    else
        error("algorithm $(algorithm) unknown")
    end
end

"""
    use_precise_ρ(cap::Intersection{N}) where {N<:Real}

Determine whether a precise algorithm for computing ``ρ`` shall be applied.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` if a precise algorithm shall be applied.

### Notes

The default implementation always returns `true`.

If the result is `false`, a coarse approximation of the support function is
returned.

This function can be overwritten by the user to control the policy.
"""
function use_precise_ρ(cap::Intersection{N}) where {N<:Real}
    return true
end

"""
    ρ(d::AbstractVector{N},
      cap::Intersection{N, S1, S2};
      [algorithm]::String="line_search",
      [kwargs...]) where {N<:Real,
                          S1<:LazySet{N},
                          S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}

Return the support function of the intersection of a compact set and a
half-space/hyperplane/line in a given direction.

### Input

- `d`         -- direction
- `cap`       -- lazy intersection of a compact set and a half-space/hyperplane/
                 line
- `algorithm` -- (optional, default: `"line_search"`): the algorithm to
                 calculate the support function; valid options are:

    * `"line_search"` -- solve the associated univariate optimization problem
                         using a line search method (either Brent or the
                         Golden Section method)
    * `"projection"`  -- only valid for intersection with a hyperplane;
                         evaluates the support function by reducing the problem
                         to the 2D intersection of a rank 2 linear
                         transformation of the given compact set in the plane
                         generated by the given direction `d` and the
                         hyperplane's normal vector `n`
    * `"simple"`      -- take the ``\\min`` of the support function evaluation
                         of each operand

### Output

The scalar value of the support function of the set `cap` in the given
direction.

### Notes

It is assumed that the set `cap.X` is compact.

Any additional number of arguments to the algorithm backend can be passed as
keyword arguments.

### Algorithm

The algorithms are based on solving the associated optimization problem

```math
\\min_{λ ∈ D_h} ρ(ℓ - λa, X) + λb.
```
where ``D_h = \\{ λ : λ ≥ 0 \\}`` if ``H`` is a half-space or
``D_h = \\{ λ : λ ∈ \\mathbb{R} \\}`` if ``H`` is a hyperplane.

For additional information we refer to:

- [G. Frehse, R. Ray. Flowpipe-Guard Intersection for Reachability Computations with
  Support Functions](https://www.sciencedirect.com/science/article/pii/S1474667015371809).
- [C. Le Guernic. Reachability Analysis of Hybrid Systems with Linear Continuous
  Dynamics, PhD thesis](https://tel.archives-ouvertes.fr/tel-00422569v2).
- [T. Rockafellar, R. Wets.
  Variational Analysis](https://www.springer.com/us/book/9783540627722).
"""
function ρ(d::AbstractVector{N},
           cap::Intersection{N, S1, S2};
           algorithm::String="line_search",
           kwargs...) where {N<:Real,
                             S1<:LazySet{N},
                             S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}
    return ρ_helper(d, cap, algorithm; kwargs...)
end

# symmetric method
function ρ(d::AbstractVector{N},
           cap::Intersection{N, S1, S2};
           algorithm::String="line_search",
           kwargs...) where {N<:Real,
                             S1<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                             S2<:LazySet{N}}
    return ρ_helper(d, swap(cap), algorithm; kwargs...)
end

"""
    ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...)
        where {N<:Real, S1<:LazySet{N}, S2<:AbstractPolyhedron{N}}

Return an upper bound on the support function of the intersection between a
compact set and a polyhedron along a given direction.

### Input

- `d`      -- direction
- `cap`    -- intersection of a compact set and a polyhedron
- `kwargs` -- additional arguments that are passed to the support-function
              algorithm

### Output

An upper bound of the support function of the given intersection.

### Algorithm

The idea is to solve the univariate optimization problem `ρ(di, X ∩ Hi)` for each
half-space in the set `P` and then take the minimum. This gives an overapproximation
of the exact support function.

This algorithm is inspired from [G. Frehse, R. Ray. Flowpipe-Guard Intersection
for Reachability Computations with Support
Functions](https://www.sciencedirect.com/science/article/pii/S1474667015371809).

### Notes

This method relies on the `constraints_list` of the polyhedron.
"""
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...
          ) where {N<:Real, S1<:LazySet{N}, S2<:AbstractPolyhedron{N}}
    return ρ_helper(d, cap; kwargs...)
end

# symmetric method
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...
          ) where {N<:Real, S1<:AbstractPolyhedron{N}, S2<:LazySet{N}}
    return ρ_helper(d, swap(cap); kwargs...)
end

# disambiguation
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...
          ) where {N<:Real, S1<:AbstractPolytope{N}, S2<:AbstractPolyhedron{N}}
    return ρ_helper(d, cap; kwargs...)
end

function ρ_helper(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...
                 ) where {N<:Real, S1<:LazySet{N}, S2<:AbstractPolyhedron{N}}
    if !use_precise_ρ(cap)
        use_simple_method = true
    else
        options = Dict(kwargs)
        use_simple_method = haskey(options, :algorithm) &&
                       options[:algorithm] == "simple"
    end

    if use_simple_method
        # simple algorithm
        return invoke(ρ, Tuple{typeof(d), Intersection{N}}, d, cap)
    end

    # more precise algorithm
    @assert isbounded(cap.X) "the first set in the intersection must be bounded"
    return minimum([ρ(d, cap.X ∩ Hi; kwargs...)
                   for Hi in constraints_list(cap.Y)])
end

"""
    ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...
     ) where {N<:Real, S1<:AbstractPolyhedron{N}, S2<:AbstractPolyhedron{N}}

Return an upper bound on the support function of the intersection between two
polyhedral sets.

### Input

- `d`      -- direction
- `cap`    -- intersection of two polyhedral sets
- `kwargs` -- additional arguments that are passed to the support-function
              algorithm

### Output

The support function for the given direction.

### Algorithm

We combine the constraints of the two polyhedra to a new `HPolyhedron`, for
which we then evaluate the support function.
"""
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}; kwargs...
          ) where {N<:Real, S1<:AbstractPolyhedron{N}, S2<:AbstractPolyhedron{N}}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

# disambiguation
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2};
           algorithm::String="line_search", kwargs...
          ) where {N<:Real, S1<:AbstractPolytope{N},
                   S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}
    return ρ_helper(d, cap, algorithm; kwargs...)
end

# symmetric method
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2};
           algorithm::String="line_search", kwargs...
          ) where {N<:Real, S1<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                   S2<:AbstractPolytope{N}}
    return ρ_helper(d, swap(cap), algorithm; kwargs...)
end

# disambiguation
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2};
           algorithm::String="line_search", kwargs...
          ) where {N<:Real, S1<:AbstractPolyhedron{N},
                   S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

# symmetric method
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2};
           algorithm::String="line_search", kwargs...
          ) where {N<:Real, S1<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                   S2<:AbstractPolyhedron{N}}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

# disambiguation
function ρ(d::AbstractVector{N}, cap::Intersection{N, S1, S2}
          ) where {N<:Real, S1<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                   S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

"""
    isbounded(cap::Intersection)

Determine whether an intersection of two convex sets is bounded.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` iff the intersection is bounded.

### Algorithm

We first check if any of the wrapped sets is bounded.
Otherwise, we check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(cap::Intersection)
    if isbounded(cap.X) || isbounded(cap.Y)
        return true
    end
    return isbounded_unit_dimensions(cap)
end

"""
    ∈(x::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}

Check whether a given point is contained in an intersection of two convex sets.

### Input

- `x`   -- point/vector
- `cap` -- intersection of two convex sets

### Output

`true` iff ``x ∈ cap``.
"""
function ∈(x::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}
    return (x ∈ cap.X) && (x ∈ cap.Y)
end

"""
    constraints_list(cap::Intersection{N}) where {N<:Real}

Return the list of constraints of an intersection of two (polyhedral) sets.

### Input

- `cap` -- intersection of two (polyhedral) sets

### Output

The list of constraints of the intersection.

### Notes

We assume that the underlying sets are polyhedral, i.e., offer a method
`constraints_list`.

### Algorithm

We create the polyhedron by taking the intersection of the `constraints_list`s of
the sets and remove redundant constraints.

This function ignores the boolean output from the in-place `remove_redundant_constraints!`,
which may inform the user that the constraints are infeasible. In that case, the
list of constraints at the moment when the infeasibility was detected is returned.
"""
function constraints_list(cap::Intersection{N}) where {N<:Real}
    constraints = [constraints_list(cap.X); constraints_list(cap.Y)]
    remove_redundant_constraints!(constraints)
    return constraints
end


# --- Intersection functions ---


"""
    isempty(cap::Intersection)

Return if the intersection is empty or not.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` iff the intersection is empty.

### Notes

The result will be cached, so a second query will be fast.
"""
function isempty(cap::Intersection)
    if isempty_known(cap)
        # use cached result
        return isempty(cap.cache)
    end
    # compute result
    empty_intersection = is_intersection_empty(cap.X, cap.Y)
    # update cache
    set_isempty!(cap, empty_intersection)

    return empty_intersection
end

"""
    plot_recipe(cap::Intersection{N}, [ε]::N=-one(N),
                [Nφ]::Int=PLOT_POLAR_DIRECTIONS) where {N<:Real}

Convert a lazy intersection to a pair `(x, y)` of points for plotting.

### Input

- `cap`  -- lazy intersection
- `ε`    -- (optional, default `0`) ignored, used for dispatch
- `Nφ`   -- (optional, default: `PLOT_POLAR_DIRECTIONS`) number of polar
            directions used in the template overapproximation

### Output

A pair `(x, y)` of points that can be plotted.
"""
function plot_recipe(cap::Intersection{N}, ε::N=zero(N),
                     Nφ::Int=PLOT_POLAR_DIRECTIONS) where {N<:Real}
    @assert dim(cap) <= 2 "cannot plot a $(dim(cap))-dimensional intersection"

    if isempty(cap)
        return plot_recipe(EmptySet{N}(dim(cap)), ε)
    elseif dim(cap) == 1
        return plot_recipe(convert(Interval, cap), ε)
    else
        # construct polygon approximation using polar directions
        P = overapproximate(cap, PolarDirections{N}(Nφ))
        return plot_recipe(P, ε)
    end
end

# ==========================================================
# Concrete operations that dispatch on a lazy Intersection
# ==========================================================

"""
    linear_map(M::AbstractMatrix{N}, cap::Intersection{N}) where {N}

Return the concrete linear map of a lazy intersection.

### Input

- `M`   -- matrix
- `cap` -- lazy intersection

### Output

The set obtained by applying the given linear map to the lazy intersection.

### Notes

This function relies on computing `cap` concretely (i.e. as a set representation),
and then applying the linear map.
"""
function linear_map(M::AbstractMatrix{N}, cap::Intersection{N}) where {N}
    return linear_map(M, intersection(cap.X, cap.Y))
end

# ==================================
# Algorithms for lazy intersection
# ==================================

function load_optim_intersection()
return quote

"""
    _line_search(ℓ, X, H::Union{<:HalfSpace, <:Hyperplane, <:Line}; [kwargs...])

Given a compact and convex set ``X`` and a halfspace ``H = \\{x: a^T x ≤ b \\}``
or a hyperplane ``H = \\{x: a^T x = b \\}``, calculate:

```math
\\min_{λ ∈ D_h} ρ(ℓ - λa, X) + λb.
```
where ``D_h = \\{ λ : λ ≥ 0 \\}`` if ``H`` is a half-space or
``D_h = \\{ λ : λ ∈ \\mathbb{R} \\}`` if ``H`` is a hyperplane.

### Input

- `ℓ`           -- direction
- `X`           -- set
- `H`           -- halfspace or hyperplane

### Output

The tuple `(fmin, λmin)`, where `fmin` is the minimum value of the function
``f(λ) = ρ(ℓ - λa) + λb`` over the feasible set ``λ ≥ 0``, and ``λmin`` is the
minimizer.

### Notes

This function requires the `Optim` package, and relies on the univariate
optimization interface `Optim.optimize(...)`.

Additional arguments to the `optimize` backend can be passed as keyword arguments.
The default method is `Optim.Brent()`.

### Examples

```jldoctest _line_search
julia> X = Ball1(zeros(2), 1.0);

julia> H = HalfSpace([-1.0, 0.0], -1.0); # x >= 0

julia> using Optim

julia> using LazySets: _line_search

julia> _line_search([1.0, 0.0], X, H) # uses Brent's method by default
(1.0, 999999.9849478417)
```

We can specify the upper bound in Brent's method:

```julia _line_search
julia> _line_search([1.0, 0.0], X, H, upper=1e3)
(1.0, 999.9999849478418)
```

Instead of using Brent, we use the Golden Section method:

```julia _line_search
julia> _line_search([1.0, 0.0], X, H, upper=1e3, method=GoldenSection())
(1.0, 381.9660112501051)
```
"""
function _line_search(ℓ, X, H::Union{<:HalfSpace, <:Hyperplane, <:Line};
                      kwargs...)
    options = Dict(kwargs)

    # Initialization
    a, b = H.a, H.b
    f(λ) = ρ(ℓ - λ[1] * a, X) + λ[1] * b

    if haskey(options, :lower)
        lower = pop!(options, :lower)
    else
        if H isa HalfSpace
            lower = 0.0
        elseif (H isa Hyperplane) || (H isa Line)
            lower = -1e6 # "big": TODO relate with f(λ)
        end
    end

    if haskey(options, :upper)
        upper = pop!(options, :upper)
    else
        upper = 1e6 # "big": TODO relate with f(λ)
    end

    if haskey(options, :method)
        method = pop!(options, :method)
    else
        method = Optim.Brent()
    end

    # Optimization
    sol = Optim.optimize(f, lower, upper, method=method, options...)

    # Recover results
    fmin, λmin = sol.minimum, sol.minimizer
    return (fmin, λmin)
end # _line_search
end # quote
end # load_optim

"""
    _projection(ℓ, X, H::Union{Hyperplane{N}, Line{N}};
                [lazy_linear_map]=false,
                [lazy_2d_intersection]=true,
                [algorithm_2d_intersection]=nothing,
                [kwargs...]) where {N}

Given a compact and convex set ``X`` and a hyperplane ``H = \\{x: n ⋅ x = γ \\}``,
calculate the support function of the intersection between the
rank-2 projection ``Π_{nℓ} X`` and the line ``Lγ = \\{(x, y): x = γ \\}``.

### Input

- `ℓ`                    -- direction
- `X`                    -- set
- `H`                    -- hyperplane
- `lazy_linear_map`      -- (optional, default: `false`) to perform the projection
                            lazily or concretely
- `lazy_2d_intersection` -- (optional, default: `true`) to perform the 2D
                            intersection between the projected set and the line
                            lazily or concretely
- `algorithm_2d_intersection` -- (optional, default: `nothing`) if given, fixes the
                                 support function algorithm used for the intersection
                                 in 2D; otherwise the default is implied

### Output

The support function of ``X ∩ H`` along direction ``ℓ``.

### Algorithm

This projection method is based on Prop. 8.2, page 103, [C. Le Guernic.
Reachability Analysis of Hybrid Systems with Linear Continuous Dynamics,
PhD thesis](https://tel.archives-ouvertes.fr/tel-00422569v2).

In the original algorithm, Section 8.2 of Le Guernic's thesis, the linear map
is performed concretely and the intersection is performed lazily (these are the
default options in this algorithm, but here the four combinations are available).
If the set ``X`` is a zonotope, its concrete projection is again a zonotope
(sometimes called "zonogon"). The intersection between this zonogon and the line
can be taken efficiently in a lazy way (see Section 8.2.2 of Le Guernic's thesis),
if one uses dispatch on `ρ(y_dir, Sℓ⋂Lγ; kwargs...)` given that `Sℓ` is itself
a zonotope.

### Notes

This function depends itself on the calculation of the support function of another
set in two dimensions. Obviously one doesn't want to use again `algorithm="projection"`
for this second calculation. The option `algorithm_2d_intersection` is such that,
if it is not given, the default support function algorithm is used (e.g. `"line_search"`).
You can still pass additional arguments to the `"line_search"` backend through the
`kwargs`.
"""
function _projection(ℓ, X, H::Union{Hyperplane{N}, Line{N}};
                     lazy_linear_map=false,
                     lazy_2d_intersection=true,
                     algorithm_2d_intersection=nothing,
                     kwargs...) where {N}

    n = H.a                  # normal vector to the hyperplane
    γ = H.b                  # displacement of the hyperplane
    Πnℓ = vcat(n', ℓ')       # projection map

    x_dir = [one(N), zero(N)]
    y_dir = [zero(N), one(N)]

    Lγ = Line(x_dir, γ)

    Xnℓ = lazy_linear_map ? LinearMap(Πnℓ, X) : linear_map(Πnℓ, X)
    Xnℓ⋂Lγ = lazy_2d_intersection ? Intersection(Xnℓ, Lγ) : intersection(Xnℓ, Lγ)

    if algorithm_2d_intersection == nothing
        return ρ(y_dir, Xnℓ⋂Lγ; kwargs...)
    else
        return ρ(y_dir, Xnℓ⋂Lγ, algorithm=algorithm_2d_intersection; kwargs...)
    end
end

"""
    get_constrained_lowdimset(cpa::CartesianProductArray{N, S},
                              P::AbstractPolyhedron{N}
                             ) where {N<:Real, S<:LazySet{N}}

Preprocess step for intersection between Cartesian product array and polyhedron.
Returns low-dimensional a `CartesianProductArray` in the constrained dimensions
of the original cpa,
constrained variables and variables in corresponding blocks, original block
structure of low-dimensional set and list of constrained blocks.

### Input

- `cpa` -- Cartesian product array of convex sets
- `P`   -- polyhedron

### Output

A tuple of low-dimensional set, list of constrained dimensions, original block
structure of low-dimensional set and corresponding blocks indices.
"""
function get_constrained_lowdimset(cpa::CartesianProductArray{N, S},
                                   P::AbstractPolyhedron{N}
                                  ) where {N<:Real, S<:LazySet{N}}

    if isbounded(P)
        blocks, non_empty_length = block_to_dimension_indices(cpa)
    else
        blocks, non_empty_length =
            block_to_dimension_indices(cpa, constrained_dimensions(P))
    end

    array = Vector{S}()
    sizehint!(array, non_empty_length)
    variables = Vector{Int}()
    block_structure = Vector{UnitRange{Int}}()
    sizehint!(block_structure, non_empty_length)

    last_var = 1
    for i in 1:length(blocks)
        start_index, end_index = blocks[i]
        block_end = last_var + end_index - start_index
        if start_index != -1
            push!(array, cpa.array[i])
            append!(variables, start_index : end_index)
            push!(block_structure, last_var : block_end)
            last_var = block_end + 1
        end
    end

    return CartesianProductArray(array), variables, block_structure, blocks
end
