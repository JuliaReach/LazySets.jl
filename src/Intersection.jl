import Base: isempty, ∈, ∩

export Intersection,
       IntersectionArray,
       array

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
    error("not implemented yet")
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
      cap::Intersection{N, <:LazySet, <:HalfSpace};
      algorithm::String="line_search",
      check_intersection::Bool=true,
      kwargs...) where {N<:AbstractFloat}

Return the support function of the intersection of a compact set and a half-space
in a given direction.

### Input

- `d`         -- direction
- `cap`       -- lazy intersection of a compact set and a half-space 
- `algorithm` -- (optional, default: `"line_search"`): the algorithm to calculate
                 the support function, valid options are:
                 
    * `"line_search"` -- solve the associated univariate optimization problem
                         using a line search method (either Brent or the
                         Golden Section method) 

- `check_intersection` -- (optional, default: `true`) if `true`, check if the 
                          intersection is empty before actually calculating the
                          support function

### Output

The scalar value of the support function of the set `cap` in the given direction.

### Notes

The `check_intersection` flag can be useful if you know in advance that the
intersection is non-empty.

Any additional number of arguments to the algorithm backend can be passed as
keyword arguments.

### Algorithm

The algorithms are based on solving the associated optimization problem

```math
\\min_\\{ λ ≥ 0 \\} ρ(ℓ - λa, X) + λb.
```

For additional information we refer to:

- [G. Frehse, R. Ray. Flowpipe-Guard Intersection for Reachability Computations with
  Support Functions](https://www.sciencedirect.com/science/article/pii/S1474667015371809).
- [C. Le Guernic. Reachability Analysis of Hybrid Systems with Linear Continuous
  Dynamics, PhD thesis](https://tel.archives-ouvertes.fr/tel-00422569v2).
- [T. Rockafellar, R. Wets.
  Variational Analysis](https://www.springer.com/us/book/9783540627722).
"""
function ρ(d::AbstractVector{N},
           cap::Intersection{N, <:LazySet, <:HalfSpace};
           algorithm::String="line_search",
           check_intersection::Bool=true,
           kwargs...) where {N<:AbstractFloat}
    
    X = cap.X    # compact set
    H = cap.Y    # halfspace

    # if the intersection is empty => stop
    if check_intersection
        is_intersection_empty(X, H) && return zero(N) # TODO use this convention? error?
    end
    
    if algorithm == "line_search"
        @assert isdefined(Main, :Optim) "the algorithm $algorithm needs " *
            "the package 'Optim' to be loaded"
        (s, _) = _line_search(d, X, H; kwargs...)
    else
        error("algorithm $(algorithm) unknown")
    end
    return s
end

function load_optim_intersection()
return quote

using Optim

"""
    _line_search(ℓ, X, H; kwargs...)

Given a compact and convex set ``X`` and a hyperplane ``H = \\{x: a^T x ≤ b \\}``,
calculate:

```math
\\min_\\{ λ ≥ 0 \\} ρ(ℓ - λa, X) + λb.
```

### Input

- `ℓ`      -- direction
- `X`      -- set
- `H`      -- hyperplane

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
function _line_search(ℓ, X, H; kwargs...)
    options = Dict(kwargs)
   
    # Initialization
    a, b = H.a, H.b
    f(λ) = ρ(ℓ - λ[1] * a, X) + λ[1] * b

    if haskey(options, :lower)
        lower = pop!(options, :lower)
    else
        lower = 0.0
    end

    if haskey(options, :upper)
        upper = pop!(options, :upper)
    else
        upper = 1e6 # TODO: relate with f(λ)
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

# Symmetric case
ρ(ℓ::AbstractVector{N}, cap::Intersection{N, <:HalfSpace, <:LazySet};
  algorithm::String="line_search", check_intersection::Bool=true,
  kwargs...) where {N<:AbstractFloat} = ρ(ℓ, cap.Y ∩ cap.X;
  algorithm=algorithm, check_intersection=check_intersection, kwargs...)
