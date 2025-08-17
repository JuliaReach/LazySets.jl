import Base: ∩

export Intersection,
       isempty_known,
       set_isempty!,
       swap,
       use_precise_ρ

"""
    IntersectionCache

Container for information cached by a lazy `Intersection` object.

### Fields

- `isempty` -- is the intersection empty? There are three possible states,
               encoded as `Int8` values -1, 0, 1:

    * ``-1`` - it is currently unknown whether the intersection is empty
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
    return c.isempty = isempty ? Int8(1) : Int8(0)
end

"""
    Intersection{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two sets.

### Fields

- `X`     -- set
- `Y`     -- set
- `cache` -- internal cache for avoiding recomputation; see
             [`IntersectionCache`](@ref)

### Notes

If the arguments of the lazy intersection are half-spaces, the set is simplified
to a polyhedron in constraint representation (`HPolyhedron`).

The intersection preserves convexity: if the set arguments are convex, then
their intersection is convex as well.

The convenience alias `∩` can be typed by `\\cap<tab>`.

### Examples

Create an expression ``Z`` that lazily represents the intersection of two
squares ``X`` and ``Y``:

```jldoctest lazy_intersection
julia> X, Y = BallInf([0.0, 0.0], 0.5), BallInf([1.0, 0.0], 0.75);

julia> Z = X ∩ Y;

julia> typeof(Z)
Intersection{Float64, BallInf{Float64, Vector{Float64}}, BallInf{Float64, Vector{Float64}}}

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
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([0.375, 0.0], [0.125, 0.5])
```
"""
struct Intersection{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
    cache::IntersectionCache

    # default constructor with dimension check
    function Intersection(X::LazySet{N}, Y::LazySet{N};
                          cache::IntersectionCache=IntersectionCache()) where {N}
        @assert dim(X) == dim(Y) "sets in an intersection must have the same " *
                                 "dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y, cache)
    end
end

# constructors simplifying to HPolyhedron
Intersection(H1::HalfSpace, H2::HalfSpace) = HPolyhedron([H1, H2])
Intersection(H::HalfSpace, P::HPolyhedron) = HPolyhedron(vcat(P.constraints, H))
Intersection(P::HPolyhedron, H::HalfSpace) = HPolyhedron(vcat(P.constraints, H))

∩(X::LazySet, Y::LazySet) = Intersection(X, Y)

isoperationtype(::Type{<:Intersection}) = true
concrete_function(::Type{<:Intersection}) = intersection

isconvextype(::Type{Intersection{N,S1,S2}}) where {N,S1,S2} = isconvextype(S1) && isconvextype(S2)

ispolyhedral(cap::Intersection) = ispolyhedral(cap.X) && ispolyhedral(cap.Y)

# Universe is the neutral element for Intersection
@neutral(Intersection, Universe)

# EmptySet is the absorbing element for Intersection
@absorbing(Intersection, EmptySet)

# interface for binary set operations
Base.first(cap::Intersection) = cap.X
second(cap::Intersection) = cap.Y
@declare_binary_operation(Intersection)

"""
    isempty_known(cap::Intersection)

Ask whether the status of emptiness is known.

### Input

- `cap` -- intersection of two sets

### Output

`true` iff the emptiness status is known.
In this case, `isempty(cap)` can be used to obtain the status in constant time.
"""
function isempty_known(cap::Intersection)
    return isempty_known(cap.cache)
end

"""
    set_isempty!(cap::Intersection, isempty::Bool)

Set the status of emptiness in the cache.

### Input

- `cap`     -- intersection of two sets
- `isempty` -- new status of emptiness
"""
function set_isempty!(cap::Intersection, isempty::Bool)
    return set_isempty!(cap.cache, isempty)
end

# equality test ignores the IntersectionCache
function ==(X::Intersection, Y::Intersection)
    if X.X != Y.X
        return false
    end
    if X.Y != Y.Y
        return false
    end
    return true
end

"""
    swap(cap::Intersection)

Return a new `Intersection` object with the arguments swapped.

### Input

- `cap` -- intersection of two sets

### Output

A new `Intersection` object with the arguments swapped.
The old cache is shared between the old and new objects.

### Notes

The advantage of using this function instead of manually swapping the arguments
is that the cache is shared.
"""
function swap(cap::Intersection)
    return Intersection(cap.Y, cap.X; cache=cap.cache)
end

"""
    dim(cap::Intersection)

Return the dimension of an intersection of two sets.

### Input

- `cap` -- intersection of two sets

### Output

The ambient dimension of the intersection of two sets.
"""
function dim(cap::Intersection)
    return dim(cap.X)
end

"""
    σ(d::AbstractVector, cap::Intersection)

Return a support vector of an intersection of two sets in a given
direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two sets

### Output

A support vector in the given direction.

### Algorithm

We compute the concrete intersection, which may be expensive.
"""
@validate function σ(d::AbstractVector, cap::Intersection)
    X = concretize(cap)
    return σ(d, X)
end

"""
    ρ(d::AbstractVector, cap::Intersection)

Return an upper bound on the support function of the intersection of two sets in
a given direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two sets

### Output

An upper bound on the support function in the given direction.

### Algorithm

The support function of an intersection of ``X`` and ``Y`` is upper-bounded by
the minimum of the support-function evaluations for ``X`` and ``Y``.
"""
@validate function ρ(d::AbstractVector, cap::Intersection)
    return _ρ_min(d, cap)
end

function _ρ_min(d::AbstractVector, cap::Intersection)
    return min(ρ(d, cap.X), ρ(d, cap.Y))
end

function ρ_helper(d::AbstractVector{M},
                  cap::Intersection{N,S1,<:Union{HalfSpace{N},Hyperplane{N},Line2D{N}}},
                  algorithm::String; kwargs...) where {M,N,S1}
    if !isbounded(cap.X)
        throw(ArgumentError("the first set in the intersection must be bounded"))
    end
    X = cap.X # compact set
    H = cap.Y # half-space or hyperplane or line

    if !use_precise_ρ(cap) || algorithm == "simple"
        return _ρ_min(d, cap)
    elseif algorithm == "line_search"
        require(@__MODULE__, :Optim; fun_name="ρ",
                explanation="(algorithm $algorithm)")
        (s, _) = _line_search(d, X, H; kwargs...)
        return s
    elseif algorithm == "projection"
        @assert H isa Hyperplane || H isa Line2D "the algorithm $algorithm " *
                                                 "cannot be used with a $(typeof(H)); it only works with hyperplanes"
        return _projection(d, X, H; kwargs...)
    else
        error("algorithm $algorithm unknown")
    end
end

"""
    use_precise_ρ(cap::Intersection)

Check whether a precise algorithm for computing ``ρ`` shall be applied.

### Input

- `cap` -- intersection of two sets

### Output

`true` if a precise algorithm shall be applied.

### Notes

The default implementation always returns `true`.

If the result is `false`, a coarse approximation of the support function is
returned.

This function can be overwritten by the user to control the policy.
"""
function use_precise_ρ(cap::Intersection)
    return true
end

"""
    ρ(d::AbstractVector, cap::Intersection{N, S1, S2};
      algorithm::String="line_search", kwargs...
     ) where {N, S1<:LazySet,
                 S2<:Union{HalfSpace, Hyperplane, Line2D}}

Evaluate the support function of the intersection of a compact set and a
half-space/hyperplane/line in a given direction.

### Input

- `d`         -- direction
- `cap`       -- lazy intersection of a compact set and a half-space/hyperplane/
                 line
- `algorithm` -- (optional, default: `"line_search"`): the algorithm to
                 calculate the support function; valid options are:

    * `"line_search"` -- solve the associated univariate optimization problem
                         using a line-search method (either Brent or the
                         Golden Section method)
    * `"projection"`  -- only valid for intersection with a hyperplane/line;
                         evaluate the support function by reducing the problem
                         to the 2D intersection of a rank-2 linear
                         transformation of the given compact set in the plane
                         generated by the given direction `d` and the
                         hyperplane's normal vector `n`
    * `"simple"`      -- take the ``\\min`` of the support-function evaluation
                         of each operand

### Output

The scalar value of the support function of the set `cap` in the given
direction.

### Notes

It is assumed that the first set of the intersection (`cap.X`) is compact.

Any additional number of arguments to the algorithm backend can be passed as
keyword arguments.

### Algorithm

The algorithms are based on solving the associated optimization problem

```math
\\min_{λ ∈ D_h} ρ(ℓ - λa, X) + λb.
```
where ``D_h = \\{ λ : λ ≥ 0 \\}`` if ``H`` is a half-space or
``D_h = \\{ λ : λ ∈ ℝ \\}`` if ``H`` is a hyperplane.

For additional information we refer to [Frehse012, LeGuernic09, RockafellarW98](@citet)
"""
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:LazySet,
                                       S2<:Union{HalfSpace,Hyperplane,Line2D}}
    return ρ_helper(d, cap, algorithm; kwargs...)
end

# symmetric method
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:Union{HalfSpace,Hyperplane,Line2D},
                                       S2<:LazySet}
    return ρ_helper(d, swap(cap), algorithm; kwargs...)
end

"""
    ρ(d::AbstractVector, cap::Intersection{N, S1, S2};
      kwargs...) where {N, S1<:LazySet, S2<:AbstractPolyhedron}

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

The idea is to solve the univariate optimization problem `ρ(di, X ∩ Hi)` for
each half-space in the polyhedron and then take the minimum. This gives an
overapproximation of the exact support value.

This algorithm is inspired from [Frehse012](@citet).

### Notes

This method relies on the `constraints_list` of the polyhedron.
"""
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     kwargs...) where {N,S1<:LazySet,S2<:AbstractPolyhedron}
    return ρ_helper(d, cap; kwargs...)
end

# symmetric method
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     kwargs...) where {N,S1<:AbstractPolyhedron,S2<:LazySet}
    return ρ_helper(d, swap(cap); kwargs...)
end

# disambiguation
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     kwargs...) where {N,S1<:AbstractPolytope,S2<:AbstractPolyhedron}
    return ρ_helper(d, cap; kwargs...)
end

function ρ_helper(d::AbstractVector, cap::Intersection{N,S1,S2};
                  kwargs...) where {N,S1<:LazySet,S2<:AbstractPolyhedron}
    if !use_precise_ρ(cap)
        use_simple_method = true
    else
        options = Dict(kwargs)
        use_simple_method = haskey(options, :algorithm) &&
                            options[:algorithm] == "simple"
    end

    if use_simple_method
        # simple algorithm
        return _ρ_min(d, cap)
    end

    # more precise algorithm
    @assert isbounded(cap.X) "the first set in the intersection must be bounded"
    return minimum([ρ(d, cap.X ∩ Hi; kwargs...)
                    for Hi in constraints_list(cap.Y)])
end

"""
    ρ(d::AbstractVector, cap::Intersection{N, S1, S2}; kwargs...
     ) where {N, S1<:AbstractPolyhedron, S2<:AbstractPolyhedron}

Evaluate the support function of the intersection between two polyhedral sets.

### Input

- `d`      -- direction
- `cap`    -- intersection of two polyhedral sets
- `kwargs` -- additional arguments that are passed to the support-function
              algorithm

### Output

The evaluation of the support function in the given direction.

### Algorithm

We combine the constraints of the two polyhedra to a new `HPolyhedron`, for
which we then evaluate the support function.
"""
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     kwargs...) where {N,S1<:AbstractPolyhedron,S2<:AbstractPolyhedron}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

# disambiguation
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:AbstractPolytope,
                                       S2<:Union{HalfSpace,Hyperplane,Line2D}}
    return ρ_helper(d, cap, algorithm; kwargs...)
end

# symmetric method
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:Union{HalfSpace,Hyperplane,Line2D},
                                       S2<:AbstractPolytope}
    return ρ_helper(d, swap(cap), algorithm; kwargs...)
end

# disambiguation
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:AbstractPolyhedron,
                                       S2<:Union{HalfSpace,Hyperplane,Line2D}}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

# symmetric method
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:Union{HalfSpace,Hyperplane,Line2D},
                                       S2<:AbstractPolyhedron}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

# disambiguation
@validate function ρ(d::AbstractVector, cap::Intersection{N,S1,S2};
                     algorithm::String="line_search",
                     kwargs...) where {N,S1<:Union{HalfSpace,Hyperplane,Line2D},
                                       S2<:Union{HalfSpace,Hyperplane,Line2D}}
    return ρ(d, HPolyhedron([constraints_list(cap.X); constraints_list(cap.Y)]))
end

"""
    isbounded(cap::Intersection)

Check whether an intersection of two sets is bounded.

### Input

- `cap` -- intersection of two sets

### Output

`true` iff the intersection is bounded.

### Algorithm

We first check if any of the wrapped sets is bounded.
Otherwise we check boundedness via
[`LazySets._isbounded_unit_dimensions`](@ref).
"""
function isbounded(cap::Intersection)
    if isbounded(cap.X) || isbounded(cap.Y)
        return true
    end
    return _isbounded_unit_dimensions(cap)
end

function isboundedtype(::Type{<:Intersection{N,S1,S2}}) where {N,S1,S2}
    return isboundedtype(S1) || isboundedtype(S2)
end

"""
    ∈(x::AbstractVector, cap::Intersection)

Check whether a given point is contained in the intersection of two sets.

### Input

- `x`   -- point/vector
- `cap` -- intersection of two sets

### Output

`true` iff ``x ∈ cap``.

### Algorithm

A point ``x`` is in the intersection iff it is in each set.
"""
@validate function ∈(x::AbstractVector, cap::Intersection)
    return (x ∈ cap.X) && (x ∈ cap.Y)
end

"""
    constraints_list(cap::Intersection)

Return a list of constraints of an intersection of two (polyhedral) sets.

### Input

- `cap` -- intersection of two (polyhedral) sets

### Output

A list of constraints of the intersection.

### Notes

We assume that the underlying sets are polyhedral, i.e., offer a method
`constraints_list`.

### Algorithm

We create the polyhedron by taking the intersection of the `constraints_list`s
of the sets and remove redundant constraints.

This function ignores the boolean output from the in-place
`remove_redundant_constraints!`, which may inform the user that the constraints
are infeasible. In that case, the list of constraints at the moment when the
infeasibility was detected is returned.
"""
function constraints_list(cap::Intersection)
    constraints = [constraints_list(cap.X); constraints_list(cap.Y)]
    remove_redundant_constraints!(constraints)
    return constraints
end

"""
    vertices_list(cap::Intersection)

Return a list of vertices of a lazy intersection of two (polyhedral) sets.

### Input

- `cap` -- intersection of two (polyhedral) sets

### Output

A list containing the vertices of the lazy intersection of two sets.

### Notes

We assume that the underlying sets are polyhedral and that the intersection is
bounded.

### Algorithm

We compute the concrete intersection using `intersection` and then take the
vertices of that representation.
"""
function vertices_list(cap::Intersection)
    return vertices_list(intersection(cap.X, cap.Y))
end

"""
    isempty(cap::Intersection)

Check whether the intersection of two sets is empty.

### Input

- `cap` -- intersection of two sets

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
    empty_intersection = isdisjoint(cap.X, cap.Y)
    # update cache
    set_isempty!(cap, empty_intersection)

    return empty_intersection
end

"""
    plot_recipe(cap::Intersection{N}, [ε]::N=-one(N),
                [Nφ]::Int=PLOT_POLAR_DIRECTIONS) where {N}

Convert an intersection of two sets to a pair `(x, y)` of points for plotting.

### Input

- `cap` -- intersection of two sets
- `ε`   -- (optional, default `0`) ignored, used for dispatch
- `Nφ`  -- (optional, default: `PLOT_POLAR_DIRECTIONS`) number of polar
           directions used in the template overapproximation

### Output

A pair `(x, y)` of points that can be plotted.
"""
function plot_recipe(cap::Intersection{N}, ε::N=zero(N),
                     Nφ::Int=PLOT_POLAR_DIRECTIONS) where {N}
    @assert dim(cap) <= 2 "cannot plot a $(dim(cap))-dimensional intersection"

    if isempty(cap)
        return plot_recipe(EmptySet{N}(dim(cap)), ε)
    elseif dim(cap) == 1
        if !isconvextype(typeof(cap))
            throw(ArgumentError("cannot plot a one-dimensional $(typeof(cap))"))
        end
        return plot_recipe(convert(Interval, cap), ε)
    else
        # construct polygon approximation using polar directions
        P = overapproximate(cap, PolarDirections{N}(Nφ))
        return plot_recipe(P, ε)
    end
end

"""
    linear_map(M::AbstractMatrix, cap::Intersection)

Return the concrete linear map of an intersection of two sets.

### Input

- `M`   -- matrix
- `cap` -- intersection of two sets

### Output

The set obtained by applying the given linear map to the intersection.

### Algorithm

This method computes the concrete intersection.
"""
@validate function linear_map(M::AbstractMatrix, cap::Intersection)
    return linear_map(M, intersection(cap.X, cap.Y))
end

"""
    _line_search(ℓ, X, H::Union{<:HalfSpace, <:Hyperplane, <:Line2D}; [kwargs...])

Given a convex set ``X`` and a half-space ``H = \\{x: a^T x ≤ b \\}`` or a
hyperplane/line ``H = \\{x: a^T x = b \\}``, calculate:

```math
\\min_{λ ∈ D_h} ρ(ℓ - λa, X) + λb.
```
where ``D_h = \\{ λ : λ ≥ 0 \\}`` if ``H`` is a half-space or
``D_h = \\{ λ : λ ∈ ℝ \\}`` if ``H`` is a hyperplane.

### Input

- `ℓ` -- direction
- `X` -- convex set
- `H` -- half-space or hyperplane or line

### Output

The tuple `(fmin, λmin)`, where `fmin` is the minimum value of the function
``f(λ) = ρ(ℓ - λa) + λb`` over the feasible set ``λ ≥ 0``, and ``λmin`` is the
minimizer.

### Notes

This function requires the `Optim` package, and relies on the univariate
optimization interface `Optim.optimize(...)`.

Additional arguments to the `optimize` backend can be passed as keyword
arguments. The default method is `Optim.Brent()`.

### Examples

```jldoctest _line_search
julia> X = Ball1(zeros(2), 1.0);

julia> H = HalfSpace([-1.0, 0.0], -1.0);  # x >= 1

julia> using Optim

julia> using LazySets: _line_search

julia> v = _line_search([1.0, 0.0], X, H);  # uses Brent's method by default

julia> v[1]
1.0
```

We can specify the upper bound in Brent's method:

```jldoctest _line_search
julia> v = _line_search([1.0, 0.0], X, H, upper=1e3);

julia> v[1]
1.0
```

Instead of Brent's method we can use the Golden Section method:

```jldoctest _line_search
julia> v = _line_search([1.0, 0.0], X, H, upper=1e3, method=GoldenSection());

julia> v[1]
1.0
```
"""
function _line_search(ℓ, X::LazySet, H::Union{<:HalfSpace,<:Hyperplane,<:Line2D};
                      kwargs...)
    return _line_search_optim(ℓ, X, H; kwargs...)
end

function load_optim_intersection()
    return quote
        function _line_search_optim(ℓ, X::LazySet, H::Union{<:HalfSpace,<:Hyperplane,<:Line2D};
                                    kwargs...)
            if !isconvextype(typeof(X))
                throw(ArgumentError("the first set in the intersection must be convex"))
            end

            options = Dict(kwargs)

            # Initialization
            a, b = H.a, H.b
            m = -ρ(-ℓ, X)  # `m` is a known lower bound for `f` below (Lemma 3 in paper)
            f(λ) = max(ρ(ℓ - λ[1] * a, X) + λ[1] * b, m)

            if haskey(options, :lower)
                lower = pop!(options, :lower)
            else
                if H isa HalfSpace
                    lower = 0.0
                elseif (H isa Hyperplane) || (H isa Line2D)
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
            sol = Optim.optimize(f, lower, upper; method=method, options...)

            # Recover results
            fmin, λmin = sol.minimum, sol.minimizer
            return (fmin, λmin)
        end
    end
end  # quote / load_optim_intersection

"""
    _projection(ℓ, X::LazySet, H::Union{Hyperplane, Line2D};
                [lazy_linear_map]=false,
                [lazy_2d_intersection]=true,
                [algorithm_2d_intersection]=nothing,
                [kwargs...])

Given a convex set ``X`` and a hyperplane ``H = \\{x: n ⋅ x = γ \\}``, calculate
the support function of the intersection between the rank-2 projection
``Π_{nℓ} X`` and the line ``Lγ = \\{(x, y): x = γ \\}``.

### Input

- `ℓ`                    -- direction
- `X`                    -- convex set
- `H`                    -- hyperplane
- `lazy_linear_map`      -- (optional, default: `false`) flag to perform the
                            projection lazily or concretely
- `lazy_2d_intersection` -- (optional, default: `true`) flag to perform the 2D
                            intersection between the projected set and the line
                            lazily or concretely
- `algorithm_2d_intersection` -- (optional, default: `nothing`) if given, fixes
                                 the support-function algorithm used for the
                                 intersection in 2D; otherwise the default is
                                 used

### Output

The evaluation of the support function of ``X ∩ H`` along direction ``ℓ``.

### Algorithm

This projection method is based on Prop. 8.2, [1, page 103].

In the original algorithm, [LeGuernic09; Section 8.2](@citet), the linear map is performed
concretely and the intersection is performed lazily (these are the default
options in this algorithm, but here the four combinations are available).
If the set ``X`` is a zonotope, its concrete projection is again a zonotope
(sometimes called "zonogon"). The intersection between this zonogon and the line
can be taken efficiently in a lazy way (see [LeGuernic09; Section 8.2.2](@cite)),
if one uses dispatch on `ρ(y_dir, Sℓ⋂Lγ; kwargs...)` given that `Sℓ` is itself
a zonotope.

### Notes

This function depends on the calculation of the support function of another set
in two dimensions. Obviously one does not want to use `algorithm="projection"`
again for this second calculation. The option `algorithm_2d_intersection` is
used for that: if not given, the default support-function algorithm is used
(e.g., `"line_search"`). You can still pass additional arguments to the
`"line_search"` backend through the `kwargs` arguments.
"""
function _projection(ℓ, X::LazySet, H::Union{Hyperplane,Line2D};
                     lazy_linear_map=false,
                     lazy_2d_intersection=true,
                     algorithm_2d_intersection=nothing,
                     kwargs...)
    if !isconvextype(typeof(X))
        throw(ArgumentError("the first set in the intersection must be convex"))
    end

    N = promote_type(eltype(X), eltype(H))
    n = H.a                  # normal vector to the hyperplane
    γ = H.b                  # displacement of the hyperplane
    Πnℓ = vcat(n', ℓ')       # projection map

    x_dir = [one(N), zero(N)]
    y_dir = [zero(N), one(N)]

    Lγ = Line2D(x_dir, γ)

    Xnℓ = lazy_linear_map ? LinearMap(Πnℓ, X) : linear_map(Πnℓ, X)
    Xnℓ⋂Lγ = lazy_2d_intersection ? Intersection(Xnℓ, Lγ) : intersection(Xnℓ, Lγ)

    if isnothing(algorithm_2d_intersection)
        return ρ(y_dir, Xnℓ⋂Lγ; kwargs...)
    else
        return ρ(y_dir, Xnℓ⋂Lγ; algorithm=algorithm_2d_intersection, kwargs...)
    end
end

"""
    get_constrained_lowdimset(cpa::CartesianProductArray{N, S},
                              P::AbstractPolyhedron{N}) where {N, S}

Preprocessing step for the intersection between a Cartesian product of a finite
number of sets and a polyhedron.

### Input

- `cpa` -- Cartesian product of a finite number of sets
- `P`   -- polyhedron

### Output

A four-tuple of:
1. a low-dimensional `CartesianProductArray` in the constrained dimensions of
   the original set `cpa`
2. the variables in the constrained blocks,
3. the original block structure of the low-dimensional sets,
4. the list of the constrained blocks.
"""
function get_constrained_lowdimset(cpa::CartesianProductArray{N,S},
                                   P::AbstractPolyhedron{N}) where {N,S}
    if isbounded(P)
        blocks, non_empty_length = block_to_dimension_indices(cpa)
    else
        blocks, non_empty_length = block_to_dimension_indices(cpa, constrained_dimensions(P))
    end

    array = Vector{S}()
    sizehint!(array, non_empty_length)
    variables = Vector{Int}()
    block_structure = Vector{UnitRange{Int}}()
    sizehint!(block_structure, non_empty_length)

    last_var = 1
    for i in eachindex(blocks)
        start_index, end_index = blocks[i]
        block_end = last_var + end_index - start_index
        if start_index != -1
            push!(array, cpa.array[i])
            append!(variables, start_index:end_index)
            push!(block_structure, last_var:block_end)
            last_var = block_end + 1
        end
    end

    return CartesianProductArray(array), variables, block_structure, blocks
end

function volume(cap::Intersection)
    return volume(intersection(cap.X, cap.Y))
end

@validate function translate(cap::Intersection, x::AbstractVector)
    X = translate(first(cap), x)
    Y = translate(second(cap), x)
    return Intersection(X, Y; cache=cap.cache)
end
