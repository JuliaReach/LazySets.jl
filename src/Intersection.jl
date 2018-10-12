import Base: isempty, ∈, ∩

export Intersection,
       IntersectionArray,
       array

global ρ_options = (false, 2, false)

"""
    Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set

### Examples

Create an expression, ``Z``, that lazily represents the intersection of two squares
``X`` and ``Y``:

```jldoctest lazy_intersection
julia> X, Y = BallInf([0,0.], 0.5), BallInf([1,0.], 0.65);

julia> Z = X ∩ Y;

julia> typeof(Z)
Intersection{Float64,BallInf{Float64},BallInf{Float64}}

julia> dim(Z)
2
```

We can check if the intersection is empty with `isempty`:

```jldoctest lazy_intersection
julia> isempty(Z)
false
```

Do not confuse `Intersection` with the concrete operation, that is computed with
the lowercase `intersection`:

```jldoctest lazy_intersection
julia> W = intersection(X, Y)
Hyperrectangle{Float64}([0.425, 0.0], [0.075, 0.5])
```
"""
struct Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function Intersection{N, S1, S2}(X::S1, Y::S2) where
            {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in an intersection must have the same " *
            "dimension"
        return new{N, S1, S2}(X, Y)
    end
end

# convenience constructor without type parameter
Intersection(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} =
    Intersection{N, S1, S2}(X, Y)

# EmptySet is the absorbing element for Intersection
@absorbing(Intersection, EmptySet)

"""
    ∩

Alias for `Intersection`.
"""
∩(X::LazySet, Y::LazySet) = Intersection(X, Y)


# --- LazySet interface functions ---


"""
    dim(cap::Intersection)::Int

Return the dimension of an intersection of two convex sets.

### Input

- `cap` -- intersection of two convex sets

### Output

The ambient dimension of the intersection of two convex sets.
"""
function dim(cap::Intersection)::Int
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
    # TODO document behavior if the direction has norm zero
    # TODO error message if the intersection is empty!
    # TODO implement
    error("the exact support vector of an intersection is not implemented yet")
end

"""
    ρ(d::AbstractVector{N}, cap::Intersection{N};
      upper_bound=false, kwargs...) where {N<:Real}

Return the support function of the intersection of two convex sets in a given
direction.

### Input

- `d`           -- direction
- `cap`         -- intersection of two convex sets
- `upper_bound` -- (optional, default: `false`) if `false`, compute the support
                   function exactly; otherwise use an overapproximative algorithm

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector{N}, cap::Intersection{N};
           upper_bound=false, kwargs...) where {N<:Real}
    if upper_bound
        return LazySets.Approximations.ρ_upper_bound(d, cap; kwargs...)
    else
        error("the exact support function of an intersection is not implemented; " *
              "try using `upper_bound=true`")
    end
end

"""
    ∈(x::AbstractVector{N}, cap::Intersection{N})::Bool where {N<:Real}

Check whether a given point is contained in an intersection of two convex sets.

### Input

- `x`   -- point/vector
- `cap` -- intersection of two convex sets

### Output

`true` iff ``x ∈ cap``.
"""
function ∈(x::AbstractVector{N}, cap::Intersection{N})::Bool where {N<:Real}
    return (x ∈ cap.X) && (x ∈ cap.Y)
end


# --- Intersection functions ---


"""
    isempty(cap::Intersection)::Bool

Return if the intersection is empty or not.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` iff the intersection is empty.
"""
function isempty(cap::Intersection)::Bool
    return is_intersection_empty(cap.X, cap.Y)
end


# ================================
# intersection of an array of sets
# ================================

"""
    IntersectionArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of a finite number of convex sets.

### Fields

- `array` -- array of convex sets

### Notes

This type assumes that the dimensions of all elements match.

The `EmptySet` is the absorbing element for `IntersectionArray`.

Constructors:

- `IntersectionArray(array::Vector{<:LazySet})` -- default constructor

- `IntersectionArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty sum with optional size hint and numeric type
"""
struct IntersectionArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

@static if VERSION < v"0.7-"
    # convenience constructor without type parameter
    IntersectionArray(arr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
        IntersectionArray{N, S}(arr)
end

# constructor for an empty sum with optional size hint and numeric type
function IntersectionArray(n::Int=0, N::Type=Float64)::IntersectionArray
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return IntersectionArray(arr)
end

# EmptySet is the absorbing element for IntersectionArray
@absorbing(IntersectionArray, EmptySet)

# add functions connecting Intersection and IntersectionArray
@declare_array_version(Intersection, IntersectionArray)

"""
    array(ia::IntersectionArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of an intersection of a finite number of convex sets.

### Input

- `ia` -- intersection of a finite number of convex sets

### Output

The array of an intersection of a finite number of convex sets.
"""
function array(ia::IntersectionArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return ia.array
end


# --- LazySet interface functions ---


"""
    dim(ia::IntersectionArray)::Int

Return the dimension of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of convex sets

### Output

The ambient dimension of the intersection of a finite number of sets.
"""
function dim(ia::IntersectionArray)::Int
    return length(ia.array) == 0 ? 0 : dim(ia.array[1])
end

"""
    σ(d::AbstractVector{<:Real}, ia::IntersectionArray)::Vector{<:Real}

Return the support vector of an intersection of a finite number of sets in a
given direction.

### Input

- `d`  -- direction
- `ia` -- intersection of a finite number of convex sets

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the individual sets.
"""
function σ(d::AbstractVector{<:Real}, ia::IntersectionArray)::Vector{<:Real}
    # TODO implement
    error("not implemented yet")
end

"""
    ∈(x::AbstractVector{N}, ia::IntersectionArray{N})::Bool where {N<:Real}

Check whether a given point is contained in an intersection of a finite number
of convex sets.

### Input

- `x`  -- point/vector
- `ia` -- intersection of a finite number of convex sets

### Output

`true` iff ``x ∈ ia``.
"""
function ∈(x::AbstractVector{N}, ia::IntersectionArray{N})::Bool where {N<:Real}
    for S in ia.array
        if x ∉ S
            return false
        end
    end
    return true
end

"""
    ρ(d::AbstractVector{N},
      cap::Intersection{N,
                        <:LazySet{N},
                        <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}};
      [algorithm]::String="line_search",
      [check_intersection]::Bool=true,
      [upper_bound]::Bool=false,
      [depth]::Int=1,
      [kwargs...]) where N<:Real

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

- `check_intersection` -- (optional, default: `true`) if `true`, check if the
                          intersection is empty before actually calculating the
                          support function
- `upper_bound` -- (optional, default: `false`) if `false`, compute the support
                   function exactly; otherwise use an overapproximative
                   algorithm
- `depth`       -- (optional, default: 1) recursion depth

### Output

The scalar value of the support function of the set `cap` in the given
direction.

### Notes

It is assumed that the set `cap.X` is compact.

The `check_intersection` flag can be useful if it is known in advance that the
intersection is non-empty.

Any additional number of arguments to the algorithm backend can be passed as
keyword arguments.

### Algorithm

The algorithms are based on solving the associated optimization problem

```math
\\min_\\{ λ ∈ D_h \\} ρ(ℓ - λa, X) + λb.
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
           cap::Intersection{N,
                             <:LazySet{N},
                             <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}};
           algorithm::String="line_search",
           check_intersection::Bool=true,
           upper_bound::Bool=false,
           depth::Int=2,
           kwargs...) where N<:Real

    X = cap.X    # compact set
    H = cap.Y    # halfspace or hyperplane

    global ρ_options
    first_execution = !ρ_options[1]
    if !first_execution
        check_intersection = false
        depth = min(ρ_options[2] - 1, -1)
        upper_bound = ρ_options[3]
    else
        depth = min(depth - 1, -1)
    end
    ρ_options = (true, depth, upper_bound)

    # if the intersection is empty => stop
    if check_intersection &&
            is_intersection_empty(X, H, false; upper_bound=upper_bound)
        error("the intersection is empty")
    end

    if upper_bound && depth < 0
        # return a naive overapproximation
        return ρ(d, cap.X)
    elseif algorithm == "line_search"
        @assert isdefined(Main, :Optim) "the algorithm $algorithm needs " *
                                        "the package 'Optim' to be loaded"
        (s, _) = _line_search(d, X, H; upper_bound=upper_bound, kwargs...)
    elseif algorithm == "projection"
        @assert H isa Hyperplane "the algorithm $algorithm cannot be used with a
                                  $(typeof(H)); it only works with hyperplanes"
        s = _projection(d, X, H; upper_bound=upper_bound, kwargs...)
    else
        error("algorithm $(algorithm) unknown")
    end
    if first_execution
        ρ_options = (false, 1, false)
    else
        ρ_options = (true, depth + 1, upper_bound)
    end
    return s
end

# Symmetric case
ρ(ℓ::AbstractVector{N},
  cap::Intersection{N,
                    <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                    <:LazySet{N}};
  algorithm::String="line_search", check_intersection::Bool=true,
  upper_bound::Bool=false, depth::Int=1, kwargs...) where N<:Real =
    ρ(ℓ, cap.Y ∩ cap.X; algorithm=algorithm,
      check_intersection=check_intersection, upper_bound=upper_bound,
      depth=depth, kwargs...)

function load_optim_intersection()
return quote

import Optim

"""
    _line_search(ℓ, X, H; [upper_bound]=false, [kwargs...])

Given a compact and convex set ``X`` and a halfspace ``H = \\{x: a^T x ≤ b \\}``
or a hyperplane ``H = \\{x: a^T x = b \\}``, calculate:

```math
\\min_\\{ λ ∈ D_h \\} ρ(ℓ - λa, X) + λb.
```
where ``D_h = \\{ λ : λ ≥ 0 \\}`` if ``H`` is a half-space or
``D_h = \\{ λ : λ ∈ \\mathbb{R} \\}`` if ``H`` is a hyperplane.

### Input

- `ℓ`           -- direction
- `X`           -- set
- `H`           -- halfspace or hyperplane
- `upper_bound` -- (optional, default: `false`) if `false`, compute the support
                   function exactly; otherwise use an overapproximative
                   algorithm
- `kwargs`      -- additional keyword arguments

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

julia> import LazySets._line_search

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
function _line_search(ℓ, X, H::Union{HalfSpace, Hyperplane, Line};
                      upper_bound::Bool=false, kwargs...)
    options = Dict(kwargs)
    # Initialization
    a, b = H.a, H.b
    ρ_rec = upper_bound ? LazySets.Approximations.ρ_upper_bound : ρ
    f(λ) = ρ_rec(ℓ - λ[1] * a, X) + λ[1] * b

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
                [upper_bound]::Bool=false,
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
- `upper_bound`          -- (optional, default: `false`) if `false`, compute the
                            support function exactly; otherwise use an
                            overapproximative algorithm

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
set in two dimensions. Obviously one doesn't want to use again `algorithm="projection"`
for this second calculation. The option `algorithm_2d_intersection` is such that,
if it is not given, the default support function algorithm is used (e.g. `"line_search"`).
You can still pass additional arguments to the `"line_search"` backend through the
`kwargs`.
"""
function _projection(ℓ, X, H::Union{Hyperplane{N}, Line{N}};
                     lazy_linear_map=false,
                     lazy_2d_intersection=true,
                     algorithm_2d_intersection=nothing,
                     upper_bound::Bool=false,
                     kwargs...) where {N}

    n = H.a                  # normal vector to the hyperplane
    γ = H.b                  # displacement of the hyperplane
    Πnℓ = vcat(n', ℓ')       # projection map

    x_dir = [one(N), zero(N)]
    y_dir = [zero(N), one(N)]

    Lγ = Line(x_dir, γ)

    Xnℓ = lazy_linear_map ? LinearMap(Πnℓ, X) : linear_map(Πnℓ, X)
    Xnℓ⋂Lγ = lazy_2d_intersection ? Intersection(Xnℓ, Lγ) : intersection(Xnℓ, Lγ)

    ρ_rec = upper_bound ? LazySets.Approximations.ρ_upper_bound : ρ
    if algorithm_2d_intersection == nothing
        return ρ_rec(y_dir, Xnℓ⋂Lγ; kwargs...)
    else
        return ρ_rec(y_dir, Xnℓ⋂Lγ, algorithm=algorithm_2d_intersection; kwargs...)
    end
end
