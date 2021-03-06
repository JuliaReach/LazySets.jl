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
    overapproximate(S::LazySet)

Alias for `overapproximate(S, Hyperrectangle)` resp. `box_approximation(S)`.
"""
overapproximate(S::LazySet) = box_approximation(S)

"""
    overapproximate(S::LazySet, ::Type{<:Hyperrectangle})

Alias for `box_approximation(S)`.
"""
function overapproximate(S::LazySet, ::Type{<:Hyperrectangle})
    return box_approximation(S)
end

# alias while Rectification is not a LazySet (#1895)
function overapproximate(r::Rectification, ::Type{<:Hyperrectangle})
    return box_approximation(r)
end

"""
    overapproximate(S::LazySet, ::Type{<:BallInf})

Alias for `ballinf_approximation(S)`.
"""
function overapproximate(S::LazySet, ::Type{<:BallInf})
    return ballinf_approximation(S)
end

"""
    overapproximate(S::LazySet{N},
                    ::Type{<:HPolygon},
                    [ε]::Real=Inf) where {N}

Return an approximation of a given 2D set using iterative refinement.

### Input

- `S`        -- convex set, assumed to be two-dimensional
- `HPolygon` -- type for dispatch
- `ε`        -- (optional, default: `Inf`) error tolerance
- `prune`    -- (optional, default: `true`) flag for removing redundant
                constraints in the end

### Output

A polygon in constraint representation.

### Notes

The result is always a convex overapproximation of the input set.

If no error tolerance ε is given, or is `Inf`, the result is a box-shaped polygon.
For convex input sets, the result is an ε-close approximation as a polygon,
with respect to the Hausdorff distance.
"""
function overapproximate(S::LazySet{N},
                         ::Type{<:HPolygon},
                         ε::Real=Inf;
                         prune::Bool=true) where {N}
    @assert dim(S) == 2 "epsilon-close approximation is only available for " *
                        "two-dimensional sets"
    if ε == Inf
        constraints = Vector{LinearConstraint{N, Vector{N}}}(undef, 4)
        constraints[1] = LinearConstraint(DIR_EAST(N), ρ(DIR_EAST(N), S))
        constraints[2] = LinearConstraint(DIR_NORTH(N), ρ(DIR_NORTH(N), S))
        constraints[3] = LinearConstraint(DIR_WEST(N), ρ(DIR_WEST(N), S))
        constraints[4] = LinearConstraint(DIR_SOUTH(N), ρ(DIR_SOUTH(N), S))
        return HPolygon(constraints, sort_constraints=false)
    else
        P = tohrep(_approximate(S, ε))
        if prune
            remove_redundant_constraints!(P)
        end
        return P
    end
end

# no-go functions for dispatch in HPolytope
function overapproximate(X::S, ::Type{<:HPolytope}) where {S<:LazySet}
    if dim(X) == 2
        throw(ArgumentError("ε-close approximation is only available using " *
              "polygons in constraint representation; try `overapproximate(X, HPolygon)`"))
    else
        throw(ArgumentError("ε-close approximation is only available for " *
              "two-dimensional sets; try `overapproximate(X, HPolytope, dirs)` " *
              "where `dirs` are template directions, e.g., `BoxDirections` " *
              "or `OctDirections`"))
    end
end

"""
    overapproximate(S::LazySet, ε::Real)

Alias for `overapproximate(S, HPolygon, ε)`.
"""
function overapproximate(S::LazySet, ε::Real)
    return overapproximate(S, HPolygon, ε)
end

# special case: overapproximation of empty set
overapproximate(∅::EmptySet, options...) = ∅

# disambiguation
overapproximate(∅::EmptySet) = ∅
for ST in LazySets.subtypes(LazySet, true)
    if ST == HPolygon  # must be defined separately below with extra argument
        continue
    end
    @eval overapproximate(∅::EmptySet, ::Type{<:$ST}) = ∅
end
overapproximate(∅::EmptySet, ::Real) = ∅
overapproximate(∅::EmptySet, ::Type{<:HPolygon}, ε::Real=Inf) = ∅
overapproximate(∅::EmptySet, ::Type{<:EmptySet}, args...) = ∅

"""
    overapproximate(X::ConvexHull{N, <:AbstractZonotope, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the convex hull of two zonotopes.

### Input

- `X`         -- convex hull of two zonotopes
- `Zonotope`  -- type for dispatch
- `algorithm` -- (optional; default: `"mean"`) choice of algorithm; possible
                 values are `"mean"` and `"join"`

### Output

A zonotope ``Z`` such that ``X ⊆ Z``.

### Algorithm

The algorithm can be controlled by the parameter `algorithm`.
Note that the results of the two implemented algorithms are generally
incomparable.

##### 'mean' method

If `algorithm == "mean"`, we choose the method proposed in [1].
The convex hull of two zonotopes ``Z₁`` and ``Z₂`` of the same order,
which we write

```math
Z_j = ⟨c^{(j)}, g^{(j)}_1, …, g^{(j)}_p⟩
```
for ``j = 1, 2``, can be overapproximated as follows:

```math
CH(Z_1, Z_2) ⊆ \\frac{1}{2}⟨c^{(1)}+c^{(2)}, g^{(1)}_1+g^{(2)}_1, …,
g^{(1)}_p+g^{(2)}_p, c^{(1)}-c^{(2)}, g^{(1)}_1-g^{(2)}_1, …, g^{(1)}_p-g^{(2)}_p⟩.
```

If the zonotope order is not the same, this algorithm calls
`reduce_order` to reduce the order to the minimum of the arguments.

It should be noted that the output zonotope is not necessarily the minimal
enclosing zonotope, which is in general expensive in high dimensions.
This is further investigated in [2].

##### 'join' method

If `algorithm == "join"`, we choose the method proposed in [3, Definition 1].
The convex hull ``X`` of two zonotopes ``Z₁`` and ``Z₂`` is overapproximated by
a zonotope ``Z₃`` such that the box approximation of ``X`` is identical with the
box approximation of ``Z₃``.
Let ``□(X)`` denote the box approximation of ``X``.
The center of ``Z₃`` is the center of ``□(X)``.

The generator construction consists of two phases.
In the first phase, we construct generators ``g`` as a combination of one
generator from ``Z₁``, say, ``g₁``, with another generator from ``Z₂``, say,
``g₂``.
The entry of ``g`` in the ``i``-th dimension is given as

```math
    g[i] = \\arg\\min_{\\min(g₁[i], g₂[i]) ≤ x ≤ \\max(g₁[i], g₂[i])} |x|.
```

If ``g`` is the zero vector, it can be omitted.

In the second phase, we construct another generator for each dimension.
These generators are scaled unit vectors.
The following formula defines the sum of all those generators.

```math
    \\sup(□(X)) - c - ∑_g |g|
```

where ``c`` is the center of the new zonotope and the ``g``s are the generators
constructed in the first phase.

##### References

[1] Reachability of Uncertain Linear Systems Using Zonotopes, A. Girard.
    HSCC 2005.

[2] Zonotopes as bounding volumes, L. J. Guibas et al, Proc. of Symposium on
    Discrete Algorithms, pp. 803-812.

[3] The zonotope abstract domain Taylor1+. K. Ghorbal, E. Goubault, S. Putot.
    CAV 2009.
"""
function overapproximate(X::ConvexHull{N, <:AbstractZonotope, <:AbstractZonotope},
                         ::Type{<:Zonotope};
                         algorithm="mean") where {N}
    # execute specific algorithm
    if algorithm == "mean"
        return _overapproximate_convex_hull_zonotope_G05(X)
    elseif algorithm == "join"
        return _overapproximate_convex_hull_zonotope_GGP09(X)
    else
        error("algorithm $algorithm is not known")
    end
end

function _overapproximate_convex_hull_zonotope_G05(X::ConvexHull{N}) where {N}
    # reduce to the same order if possible
    m1, m2 = ngens(X.X), ngens(X.Y)
    if m1 < m2
        Y = CH(X.X, reduce_order(X.Y, m1))
    elseif m1 > m2
        Y = CH(reduce_order(X.X, m2), X.Y)
    else
        Y = X
    end
    Z1, Z2 = Y.X, Y.Y

    if order(Z2) > order(Z1)
        Z1, Z2 = Z2, Z1
    end

    c1 = center(Z1)
    c2 = center(Z2)
    G1 = genmat(Z1)
    G2 = genmat(Z2)
    c = (c1 + c2) / N(2)

    # the case of equal order is treated separately to avoid a slicing
    # (this creates a copy)
    if order(Z1) == order(Z2)
        G = hcat(G1 .+ G2,
                 c1 - c2,
                 G1 .- G2) / N(2)
    else
        G = hcat((G1[:, 1:ngens(Z2)] .+ G2) / N(2),
                 (c1 - c2) / N(2),
                 (G1[:, 1:ngens(Z2)] .- G2) / N(2),
                 G1[:, ngens(Z2)+1:end])
    end
    Z = Zonotope(c, G)
    return remove_zero_generators(Z)
end

function _overapproximate_convex_hull_zonotope_GGP09(X::ConvexHull{N}) where {N}
    Z1, Z2 = X.X, X.Y
    m = min(ngens(Z1), ngens(Z2))
    G1, G2 = genmat(Z1), genmat(Z2)
    n = dim(Z1)
    box = box_approximation(X)

    # new center: mid point of box approximation
    c = center(box)

    # the k-th new generator is a simple combination of the old k-th generators
    G = Vector{Vector{N}}()
    sizehint!(G, m + n)
    g_sum = zeros(N, n)
    @inbounds for j in 1:m
        g1, g2 = G1[:, j], G2[:, j]
        g = Vector{N}(undef, n)
        for i in 1:n
            gi_min, gi_max = g1[i] < g2[i] ? (g1[i], g2[i]) : (g2[i], g1[i])
            if gi_min <= zero(N) && gi_max >= zero(N)
                g[i] = zero(N)
            elseif abs(gi_min) <= abs(gi_max)
                g[i] = gi_min
            else
                g[i] = gi_max
            end
        end
        if !iszero(g)
            push!(G, g)
            g_sum += abs.(g)
        end
    end

    # one more new generator (a scaled unit vector) for every dimension
    g_total = high(box) - c - g_sum
    for i in 1:n
        if !iszero(g_total[i])
            g = SingleEntryVector(i, n, g_total[i])
            push!(G, g)
        end
    end

    return Zonotope(c, G)
end

"""
    overapproximate(lm::LinearMap{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate a lazy linear map of a zonotopic set with a zonotope.

### Input

- `lm`       -- lazy linear map of a zonotopic set
- `Zonotope` -- type for dispatch

### Output

The tight zonotope corresponding to `lm`.
"""
function overapproximate(lm::LinearMap{N, <:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    return convert(Zonotope, lm)
end

"""
    overapproximate(X::LazySet{N}, dir::AbstractDirections; [prune]::Bool=true) where {N}

Overapproximate a (possibly unbounded) set with template directions.

### Input

- `X`     -- set
- `dir`   -- (concrete) direction representation
- `prune` -- (optional, default: `true`) flag for removing redundant
             constraints in the end

### Output

A polyhedron overapproximating the set `X` with the directions from `dir`.
The overapproximation is computed using support functions. If the obtained set is
bounded, the result is an `HPolytope`. Otherwise the result is an `HPolyhedron`.
"""
function overapproximate(X::LazySet{N}, dir::AbstractDirections{N, VN}; prune::Bool=true) where {N, VN}
    n = dim(X)
    H = Vector{LinearConstraint{N, VN}}()
    sizehint!(H, length(dir))

    for d in dir
        sf = ρ(d, X)
        if !isinf(sf)
            push!(H, LinearConstraint(d, sf))
        end
    end
    if prune
        remove_redundant_constraints!(H) || throw(ArgumentError("unable to remove redundant constraints"))
    end

    # if the input is bounded and the directions are bounding => output is bounded
    # otherwise, check boundedness of the output
    if (isbounded(X) && isbounding(dir)) || _isbounded_stiemke(HPolyhedron(H))
        return HPolytope(H, check_boundedness=false)
    else
        return HPolyhedron(H)
    end
end

# alias with HPolytope type as second argument
function overapproximate(X::LazySet{N}, ::Type{<:HPolytope}, dirs::AbstractDirections{N}; prune::Bool=true) where {N}
    P = overapproximate(X, dirs, prune=prune)
    isbounded(P) || throw(ArgumentError("can't overapproximate with an `HPolytope` " *
                                        "because the set is unbounded; try using an `HPolyhedron`"))
    return convert(HPolytope, P)
end

# alias with HPolyhedron type as second argument
function overapproximate(X::LazySet{N}, ::Type{<:HPolyhedron}, dirs::AbstractDirections{N}; prune::Bool=true) where {N}
    return convert(HPolyhedron, overapproximate(X, dirs, prune=true))
end

# disambiguation
overapproximate(∅::EmptySet{N}, ::Type{<:HPolytope}, dirs::AbstractDirections{N};
                prune::Bool=true) where {N} = ∅
overapproximate(∅::EmptySet{N}, ::Type{<:HPolyhedron},
                dirs::AbstractDirections{N}; prune::Bool=true) where {N} = ∅

# this function overapproximates a bounded polyhedron with a list of directions
# that define a bounded set (without checking these assumptions); the result is always bounded
function _overapproximate_bounded_polyhedron(X::LazySet{N}, dir::AbstractDirections{N, VN}; prune::Bool=true) where {N, VN}
    H = Vector{LinearConstraint{N, VN}}()
    sizehint!(H, length(dir))
    for d in dir
        sf = ρ(d, X)
        push!(H, LinearConstraint(d, sf))
    end
    if prune
        remove_redundant_constraints!(H) || throw(ArgumentError("unable to remove redundant constraints"))
    end
    return HPolytope(H, check_boundedness=false)
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
                         dir::Type{<:AbstractDirections}; kwargs...) where {N}
    return overapproximate(X, dir{N}(dim(X)); kwargs...)
end

# disambiguation
overapproximate(∅::EmptySet, dir::Type{<:AbstractDirections}; kwargs...) = ∅
overapproximate(∅::EmptySet, dir::AbstractDirections; prune::Bool=true) = ∅

"""
    overapproximate(S::LazySet{N}, ::Type{<:Interval}) where {N}

Return the overapproximation of a unidimensional set with an interval.

### Input

- `S`        -- one-dimensional set
- `Interval` -- type for dispatch

### Output

An interval.

### Algorithm

We use two support-function evaluations.
"""
function overapproximate(S::LazySet{N}, ::Type{<:Interval}) where {N}
    @assert dim(S) == 1 "cannot overapproximate a $(dim(S))-dimensional set " *
                        "with an `Interval`"
    return Interval(-ρ(N[-1], S), ρ(N[1], S))
end

"""
    overapproximate(cap::Intersection, ::Type{<:Interval})

Return the overapproximation of a unidimensional intersection with an interval.

### Input

- `cap`      -- one-dimensional lazy intersection
- `Interval` -- type for dispatch

### Output

An interval.

### Algorithm

The algorithm recursively overapproximates the two intersected sets with
intervals and then intersects these.
"""
function overapproximate(cap::Intersection, ::Type{<:Interval})
    @assert dim(cap) == 1 "cannot overapproximate a $(dim(cap))-dimensional " *
                        "intersection with an `Interval`"
    X = overapproximate(cap.X, Interval)
    Y = overapproximate(cap.Y, Interval)
    return intersection(X, Y)
end

"""
    overapproximate(cap::IntersectionArray, ::Type{<:Interval})

Return the overapproximation of a unidimensional intersection with an interval.

### Input

- `cap`      -- one-dimensional lazy intersection
- `Interval` -- type for dispatch

### Output

An interval.

### Algorithm

The algorithm recursively overapproximates the two intersected sets with
intervals and then intersects these.
"""
function overapproximate(cap::IntersectionArray, ::Type{<:Interval})
    @assert dim(cap) == 1 "cannot overapproximate a $(dim(cap))-dimensional " *
                        "intersection with an `Interval`"
    a = array(cap)
    X = overapproximate(a[1], Interval)
    @inbounds for i in 2:length(a)
        Y = a[i]
        X = intersection(X, Y)
        if isempty(X)
            return X
        end
    end
    return X
end

function overapproximate_cap_helper(X::LazySet,             # convex set
                                    P::AbstractPolyhedron,  # polyhedron
                                    dir::AbstractDirections;
                                    kwargs...
                                   )
    Hi = constraints_list(P)
    m = length(Hi)
    N = promote_type(eltype(X), eltype(P))
    constraints = Vector{HalfSpace{N, Vector{N}}}() # TODO: use directions type, see #2031
    sizehint!(constraints, length(dir))
    return_type = HPolytope

    for di in dir
        ρ_X_Hi_min = ρ(di, X ∩ Hi[1], kwargs...)
        for i in 2:m
            ρ_X_Hi = ρ(di, X ∩ Hi[i], kwargs...)
            if ρ_X_Hi < ρ_X_Hi_min
                ρ_X_Hi_min = ρ_X_Hi
            end
        end
        if ρ_X_Hi_min == N(Inf)
            # unbounded in this direction => return a polyhedron later
            return_type = HPolyhedron
        else
            push!(constraints, _normal_Vector(HalfSpace(di, ρ_X_Hi_min))) # TODO: remove vector
        end
    end
    return return_type(constraints)
end

"""
    overapproximate(cap::Intersection{N, <:LazySet, <:AbstractPolyhedron},
                    dir::AbstractDirections;
                    kwargs...
                   ) where {N}

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
                                           <:AbstractPolyhedron},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:LazySet},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:AbstractPolyhedron},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    # important: the result may not be a polytope!
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:AbstractPolyhedron},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:AbstractPolytope},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.Y, cap.X, dir; kwargs...)
end

# disambiguation
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:AbstractPolytope},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    return overapproximate_cap_helper(cap.X, cap.Y, dir; kwargs...)
end

"""
    overapproximate(cap::Intersection{N, <:HalfSpace, <:AbstractPolytope},
                    dir::AbstractDirections;
                    [kwargs]...
                   ) where {N}

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
                                           <:HalfSpace,
                                           <:AbstractPolytope},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
    H = HPolytope{N, Vector{N}}()
    c = H.constraints
    push!(c, _normal_Vector(cap.X))
    append!(c, _normal_Vector(cap.Y))
    return overapproximate(H, dir; kwargs...)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolytope,
                                           <:HalfSpace},
                         dir::AbstractDirections;
                         kwargs...
                        ) where {N}
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
                     constant_term, evaluate, mid, get_numvars

@inline function get_linear_coeffs(p::Taylor1)
    if p.order == 0
        return zeros(eltype(p), 1)
    end
    return linear_polynomial(p).coeffs[2:end]
end

@inline function get_linear_coeffs(p::TaylorN)
    if p.order == 0
        n = get_numvars()
        return zeros(eltype(p), n)
    end
    return linear_polynomial(p).coeffs[2].coeffs
end

"""
    overapproximate(vTM::Vector{TaylorModel1{T, S}},
                    ::Type{<:Zonotope}) where {T, S}

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
                            ::Type{<:Zonotope}) where {T, S}
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
                    ::Type{<:Zonotope}) where {N,T, S}


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
                         ::Type{<:Zonotope}) where {N, T, S}
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
    Z = Zonotope(c, hcat(gen_lin, Diagonal(gen_rem)))
    return remove_zero_generators(Z)
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
    overapproximate(lm::LinearMap{N, <:AbstractZonotope, NM,
                                  <:AbstractIntervalMatrix{NM}},
                    ::Type{<:Zonotope}) where {N, NM}

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
function overapproximate(lm::LinearMap{N, <:AbstractZonotope, NM,
                                       <:AbstractIntervalMatrix{NM}},
                         ::Type{<:Zonotope}) where {N, NM}
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
    overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                    ::Type{CartesianProductArray{N, S}}
                   ) where {N, S<:LazySet}

Decompose a lazy linear map of a Cartesian product array while keeping the
original block structure.

### Input

- `lm`                    -- lazy linear map of Cartesian product array
- `CartesianProductArray` -- type for dispatch

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                         ::Type{CartesianProductArray{N, S}}
                        ) where {N, S<:LazySet}
    cpa = array(lm.X)
    arr = Vector{S}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, S)
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray},
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
function overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                         ::Type{<:CartesianProductArray},
                         dir::Type{<:AbstractDirections}) where {N}
    cpa = array(lm.X)
    arr = Vector{HPolytope{N}}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, dir)
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray},
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
function overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                         ::Type{<:CartesianProductArray},
                         set_type::Type{<:LazySet}) where {N}
    cpa = array(lm.X)
    arr = Vector{set_type{N}}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, set_type)
end

function _overapproximate_lm_cpa!(arr, M, cpa, overapprox_option)
    # construct Minkowski sum for block row
    function _block_row(cpa::Vector{S}, M::AbstractMatrix,
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
    overapproximate(rm::ResetMap{N, <:CartesianProductArray},
                    ::Type{<:CartesianProductArray}, oa) where {N}

Overapproximate a reset map (that only resets to zero) of a Cartesian product
by a new Cartesian product.

### Input

- `rm`                    -- reset map
- `CartesianProductArray` -- type for dispatch
- `oa`                    -- overapproximation option

### Output

A Cartesian product with the same block structure.

### Notes

This implementation currently only supports resets to zero.

### Algorithm

We convert the `ResetMap` into a `LinearMap` and then call the corresponding
overapproximation method.
"""
function overapproximate(rm::ResetMap{N, <:CartesianProductArray},
                         ::Type{<:CartesianProductArray}, oa) where {N}
    if any(!iszero, values(rm.resets))
        error("we currently only support resets to zero")
    end

    lm = matrix(rm) * rm.X
    return overapproximate(lm, CartesianProductArray, oa)
end

"""
    overapproximate(cap::Intersection{N,
                                      <:CartesianProductArray,
                                      <:AbstractPolyhedron},
                    ::Type{CartesianProductArray}, oa) where {N}

Return the intersection of the Cartesian product of a finite number of convex
sets and a polyhedron.

### Input

- `cap`                   -- lazy intersection of a Cartesian product array and a polyhedron
- `CartesianProductArray` -- type for dispatch
- `oa`                    -- overapproximation option

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
                                           <:CartesianProductArray,
                                           <:AbstractPolyhedron},
                            ::Type{CartesianProductArray}, oa) where {N}

    cpa, P = cap.X, cap.Y

    cpa_low_dim, vars, block_structure, blocks =
        get_constrained_lowdimset(cpa, P)

    hpoly_low_dim = HPolytope(constraints_list(cpa_low_dim))
    low_intersection = intersection(hpoly_low_dim, project(P, vars))

    if isempty(low_intersection)
        return EmptySet{N}(dim(cap))
    end

    decomposed_low_set = decompose(low_intersection, block_structure, oa)

    return substitute_blocks(decomposed_low_set, cpa, blocks)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                            <:AbstractPolyhedron,
                                            <:CartesianProductArray},
                            ::Type{CartesianProductArray}, oa) where {N}
    overapproximate(Intersection(cap.Y, cap.X), oa)
end

"""
    overapproximate(Z::Zonotope{N}, ::Type{<:Zonotope}, r::Union{Integer, Rational}) where {N}

Reduce the order of a zonotope by overapproximating with a zonotope with less
generators.

### Input

- `Z` -- zonotope
- `Zonotope` -- desired type for dispatch
- `r` -- desired order

### Output

A new zonotope with less generators, if possible.

### Algorithm

This function implements the algorithm described in A. Girard's
*Reachability of Uncertain Linear Systems Using Zonotopes*, HSCC. Vol. 5. 2005.

If the desired order is smaller than one, the zonotope is *not* reduced.
"""
function overapproximate(Z::Zonotope{N}, ::Type{<:Zonotope}, r::Union{Integer, Rational}) where {N}
    c, G = Z.center, Z.generators
    d, p = dim(Z), ngens(Z)

    if r * d >= p || r < 1
        # do not reduce
        return Z
    end

    h = zeros(N, p)
    for i in 1:p
        h[i] = norm(G[:, i], 1) - norm(G[:, i], Inf)
    end
    ind = sortperm(h)

    m = p - floor(Int, d * (r - 1)) # subset of ngens that are reduced
    rg = G[:, ind[1:m]] # reduced generators

    # interval hull computation of reduced generators
    Gbox = Diagonal(sum(abs.(rg), dims=2)[:])
    if m < p
        Gnotred = G[:, ind[m+1:end]]
        Gred = [Gnotred Gbox]
    else
        Gred = Gbox
    end
    return Zonotope(c, Gred)
end

"""
    overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                    dir::AbstractDirections;
                    [algorithm]="vrep", kwargs...)

Overapproximate a polytopic set with a zonotope.

### Input

- `X`         -- polytopic set
- `Zonotope`  -- type for dispatch
- `dir`       -- directions used for the generators
- `algorithm` -- (optional, default: `"vrep"`) method used to compute the overapproximation
- `kwargs`    -- further algorithm choices

### Output

A zonotope that overapproximates `X` and uses at most the directions provided in
`dir` (redundant directions will be ignored).

### Notes

Two algorithms are available:

- `"vrep"` -- Overapproximate a polytopic set with a zonotope of minimal total generator sum
              using only generators in the given directions. Under this constraint,
              the zonotope has the minimal sum of generator vectors. See the docstring
              of [`_overapproximate_zonotope_vrep`](@ref) for further details.

- `"cpa"` -- Overapproximate a polytopic set with a zonotope using a cartesian
             decomposition into two-dimensional blocks. See the docstring
             of [`_overapproximate_zonotope_cpa`](@ref) for further details.
"""
function overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                         dir::AbstractDirections;
                         algorithm="vrep", kwargs...)
    if algorithm == "vrep"
        return  _overapproximate_zonotope_vrep(X, dir, kwargs...)
    elseif algorithm == "cpa"
        cpa = _overapproximate_zonotope_cpa(X, dir)
        return convert(Zonotope, cpa)
    else
        throw(ArgumentError("algorithm $algorithm is not known"))
    end
end

function overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                         dir::Type{<:AbstractDirections};
                         algorithm="vrep", kwargs...)
    overapproximate(X, ZT, dir(dim(X)), algorithm=algorithm, kwargs...)
end

# disambiguation
function overapproximate(∅::EmptySet, ZT::Type{<:Zonotope},
                         dir::AbstractDirections;
                         algorithm="vrep", kwargs...)
    return ∅
end
function overapproximate(∅::EmptySet, ZT::Type{<:Zonotope},
                         dir::Type{<:AbstractDirections};
                         algorithm="vrep", kwargs...)
    return ∅
end

"""
    _overapproximate_zonotope_vrep(X::LazySet{N},
                                   dir::AbstractDirections;
                                   solver=default_lp_solver(N)) where {N}

Overapproximate a polytopic set with a zonotope of minimal total generator sum
using only generators in the given directions.

### Input

- `X`        -- polytopic set
- `dir`      -- directions used for the generators
- `solver`   -- (optional, default: `default_lp_solver(N)`) the backend used to
                solve the linear program

### Output

A zonotope that overapproximates `X` and uses at most the directions provided in
`dir` (redundant directions will be ignored).
Under this constraint, the zonotope has the minimal sum of generator vectors.

### Notes

The algorithm only requires one representative of each generator direction and
their additive inverse (e.g. only one of `[1, 0]` and `[-1, 0]`) and assumes
that the directions are normalized.
We preprocess the directions in that respect.

### Algorithm

We solve a linear program parametric in the vertices ``v_j`` of `X` and the
directions ``d_k`` in `dir` presented in Section 4.2 in [1], adapting the
notation to the one used in this library.

```math
    \\min \\sum_{k=1}^l α_k \\
    s.t. \\
    c + \\sum_{k=1}^l b_{kj} * d_k = v_j \\quad \\forall j \\
    -α_k ≤ b_{kj} ≤ α_k \\quad \\forall k, j \\
    α_k ≥ 0 \\quad \\forall k
```

The resulting zonotope has center `c` and generators `α_k · d_k`.

Note that the first type of side constraints is vector-based and that the
nonnegativity constraints (last type) are not stated explicitly in [1].

[1] Zonotopes as bounding volumes, L. J. Guibas et al, Proc. of Symposium on
    Discrete Algorithms, pp. 803-812.
"""
function _overapproximate_zonotope_vrep(X::LazySet{N},
                                        dir::AbstractDirections;
                                        solver=default_lp_solver(N)) where {N}
    # TODO "normalization" here involves two steps: removing opposite directions
    # and normalizing the direction vector
    # for the latter we can use the normalization information from dispatch on
    # the directions (some typical such directions are normalized by definition)
    function normalize_directions(dir::AD) where {AD}
        dirs = Vector{Vector{N}}()
        for d in dir
            d_normalized = normalize(d)
            if (d_normalized ∉ dirs) && (-d_normalized ∉ dirs)
                push!(dirs, d_normalized)
            end
        end
        return dirs
    end

    n = dim(X)
    dirs = normalize_directions(dir)
    l = length(dirs)
    V = vertices_list(X)
    m = length(V)
    nvariables = l + n + l * m  # l 'α_k' + n 'p[i]' + lm 'b_lj'
    nconstraints = n * m + 2 * l * m  # nm 'p_j' + 2lm 'b_kj'

    obj = vcat(ones(l), zeros(nvariables - l))
    A = spzeros(N, nconstraints, nvariables)
    sense = vcat(fill('=', n * m), fill('<', 2 * l * m))
    b = zeros(N, nconstraints)
    r = 1
    # constraints p[i] + \sum_k v_j*b_kj = p_j
    col_offset_p = l
    col_offset_b = l + n + 1
    for (j, vj) in enumerate(V)
        for i in 1:n
            A[r, col_offset_p + i] = N(1)  # p[i]
            for (k, dk) in enumerate(dirs)
                A[r, col_offset_b + (k - 1) * m] = dk[i]
            end
            b[r] = vj[i]
            r += 1
        end
        col_offset_b += 1
    end
    # constraints -α_k ± b_kj <= 0
    col_offset_b = l + n
    for k in 1:l
        for j in 1:m
            for sign_b in N[1, -1]
                A[r, k] = N(-1)  # - α_k
                A[r, col_offset_b + j] = sign_b  # ± b_kj
                r += 1
            end
        end
        col_offset_b += m
    end
    lower_bounds = vcat(zeros(N, l), fill(N(-Inf), nvariables - l))
    lp = linprog(obj, A, sense, b, lower_bounds, Inf, solver)

    if lp.status != :Optimal
        error("got unexpected status from LP solver: $(lp.status)")
    end
    c = lp.sol[(col_offset_p+1):(col_offset_p+n)]
    ck = lp.sol[1:l]
    G = Matrix{N}(undef, n, l)
    for (j, dj) in enumerate(dirs)
        G[:, j] = ck[j] * dj
    end
    return Zonotope(c, G)
end

# overload on direction type
function _overapproximate_zonotope_vrep(X::LazySet{N},
                                        dir::Type{<:AbstractDirections};
                                        solver=default_lp_solver(N)) where {N}
    return _overapproximate_zonotope_vrep(X, dir(dim(X)), solver=solver)
end

"""
    _overapproximate_zonotope_cpa(X::LazySet, dir::Type{<:AbstractDirections})

Overapproximate a polytopic set with a zonotope using cartesian decomposition.

### Input

- `X`        -- polytopic set
- `dir`      -- directions used for the generators

### Output

A zonotope that overapproximates `X`.

### Notes

The algorithm decomposes `X` in 2D sets and overapproximates those sets with
zonotopes, and finally takes the cartesian product of the sets and converts to a zonotope.

### Algorithm

The algorithm used is based on the section 8.2.4 of [1].

[1] Le Guernic, C. (2009). Reachability analysis of hybrid systems with linear
continuous dynamics (Doctoral dissertation).
"""
function _overapproximate_zonotope_cpa(X::LazySet, dir::Type{<:AbstractDirections})
    n = dim(X)
    # overapproximate 2D blocks
    if n > 1
        πX_2D = [project(X, [i, i+1]) for i in 1:2:n-1]
        Z_2D = [_overapproximate_zonotope_vrep(poly, dir(2)) for poly in πX_2D]
    end

    if iseven(n)
        out = Z_2D
    else
        # odd case projects onto an interval
        πX_n = project(X, [n])
        πX_n = overapproximate(πX_n, Interval)
        πX_n = convert(Zonotope, πX_n)

        if n == 1
            out = [πX_n]
        else
            out = vcat(Z_2D, πX_n)
        end
    end
    return CartesianProductArray(out)
end

# overload on direction type
function _overapproximate_zonotope_cpa(X::LazySet,
                                       dir::AbstractDirections)
    _overapproximate_zonotope_cpa(X, typeof(dir))
end

"""
    overapproximate(r::Rectification{N, <:AbstractZonotope}, ::Type{<:Zonotope}) where {N}

Overapproximation of the rectification of a zonotopic set.

### Input

- `r` -- lazy rectification of a zonotopic set
- `Zonotope` -- type for dispatch

### Output

A zonotope overapproximation of the set obtained by rectifying `Z`.

### Algorithm

This function implements [Theorem 3.1, 1].

[1] *Singh, G., Gehr, T., Mirman, M., Püschel, M., & Vechev, M. (2018). Fast
and effective robustness certification. In Advances in Neural Information
Processing Systems (pp. 10802-10813).*
"""
function overapproximate(r::Rectification{N, <:AbstractZonotope}, ::Type{<:Zonotope}) where {N}
    Z = copy(set(r))
    c = center(Z)
    G = genmat(Z)
    n, m = size(G)
    H = overapproximate(Z, Hyperrectangle)
    row_idx = Vector{Int}()
    μ_idx = Vector{N}()

    @inbounds for i in 1:n
        lx, ux = low(H, i), high(H, i)
        if !_leq(lx, zero(N))
            nothing
        elseif _leq(ux, zero(N))
            c[i] = zero(N)
            for j in 1:m
                G[i, j] = zero(N)
            end
        else
            λ = ux / (ux - lx)
            μ = - λ * lx / 2
            c[i] = c[i] * λ + μ
            for j in 1:m
                G[i, j] = G[i, j] * λ
            end
            push!(row_idx, i)
            push!(μ_idx, μ)
        end
    end

    q = length(row_idx)
    if q >= 1
        Gnew = zeros(N, n, q)
        j = 1
        @inbounds for i in row_idx
            Gnew[i, j] = μ_idx[j]
            j += 1
        end
        Gout = hcat(G, Gnew)
    else
        Gout = G
    end

    return Zonotope(c, remove_zero_columns(Gout))
end

"""
    overapproximate(CHA::ConvexHullArray{N, <:AbstractZonotope}, ::Type{<:Zonotope}) where {N}

Overapproximation of the convex hull array of zonotopic sets.

### Input

- `CHA` -- convex hull array of zonotopic sets
- `Zonotope` -- type for dispatch

### Output

A zonotope overapproximation of the convex hull array of zonotopic sets.

### Algorithm

This function iteratively applies the overapproximation algorithm for the
convex hull of two zonotopes to the given array of zonotopes.
"""
function overapproximate(CHA::ConvexHullArray{N, <:AbstractZonotope}, ::Type{<:Zonotope}) where {N}
    arr = array(CHA)
    n = length(arr)
    if n == 1
        return arr[1]
    else
        Zaux = overapproximate(ConvexHull(arr[1], arr[2]), Zonotope)
        for k in 3:n
            Zaux = overapproximate(ConvexHull(Zaux, arr[k]), Zonotope)
        end
        return Zaux
    end
end

"""
    overapproximate(Z::AbstractZonotope, ::Type{<:HParallelotope}, indices=1:dim(Z))

Overapproximation of a zonotopic set with a parallelotopic set in constraint
representation.

### Input

- `Z`              -- zonotopic set
- `HParallelotope` -- type for dispatch
- `indices`        -- (optional; default: `1:dim(Z)`) generator indices selected
                       when constructing the parallelotope

### Output

An overapproximation of the given zonotope using a parallelotope.

### Algorithm

The algorithm is based on Proposition 8 discussed in Section 5 of [1].

[1] Althoff, M., Stursberg, O., & Buss, M. (2010). *Computing reachable sets of
hybrid systems using a combination of zonotopes and polytopes*. Nonlinear
analysis: hybrid systems, 4(2), 233-249.
"""
function overapproximate(Z::AbstractZonotope, ::Type{<:HParallelotope}, indices=1:dim(Z))
    Zred = _overapproximate_hparallelotope(Z, indices)
    return convert(HParallelotope, Zred)
end

function _overapproximate_hparallelotope(Z::AbstractZonotope, indices=1:dim(Z))
    length(indices) == dim(Z) || throw(ArgumentError("the number of generator indices is $(length(indices)), " *
                                                     "but it was expected to be $(dim(Z))"))

    p, n = ngens(Z), dim(Z)
    if p == n
        return Z
    elseif p < n
        error("the zonotope order is $(order(Z)) but it should be at least 1")
    end

    G = genmat(Z)
    Γ = G[:, indices]
    □Γ⁻¹Z = box_approximation(linear_map(inv(Γ), Z))
    return linear_map(Γ, □Γ⁻¹Z)
end

"""
    overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                    dirs::AbstractDirections) where {N}

Overapproximation of the intersection between a zonotopic set and a hyperplane

### Input

- `X`    -- intersection between a zonotopic set and a hyperplane
- `dirs` -- type of direction representation

### Output

An overapproximation of the intersection between a zonotopic set and a hyperplane.

### Algorithm

This function implements [Algorithm 8.1, 1].

[1] *Colas Le Guernic. Reachability Analysis of Hybrid Systems with Linear
Continuous Dynamics. Computer Science [cs]. Université Joseph-Fourier - Grenoble
I, 2009. English. fftel-00422569v2f*
"""
function overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                         dirs::AbstractDirections) where {N}
    dim(X) == dim(dirs) || throw(ArgumentError("the dimension of the set, $(dim(X)) doesn't" *
                                 " match the dimension of the template, $(dim(dirs))"))
    Z, G = X.X, X.Y

    if isdisjoint(Z, G)
        return EmptySet{N}(dim(Z))
    end

    n = G.a                           # normal vector to the hyperplane
    γ = G.b                           # displacement of the hyperplane
    Lᵧ = Line2D([one(N), zero(N)], γ)  # line (x, y) : x = γ

    constraints = Vector{HalfSpace{N, eltype(dirs)}}()
    for l in dirs
        Πₙₗ = vcat(n', l')              # projection map
        πZₙₗ = linear_map(Πₙₗ, Z)

        ρₗ = LazySets._bound_intersect_2D(πZₙₗ, Lᵧ)

        push!(constraints, HalfSpace(l, ρₗ))
    end
    T = isbounding(dirs) ? HPolytope : HPolyhedron
    return T(constraints)
end

function overapproximate(X::Intersection{N, <:Hyperplane, <:AbstractZonotope},
                         dirs::AbstractDirections) where {N}
    return overapproximate(X.Y ∩ X.X, dirs)
end

function overapproximate(X::Intersection{N, <:AbstractZonotope,
                         <:Hyperplane}, dirs::Type{<:AbstractDirections}) where {N}
    return overapproximate(X, dirs(dim(X)))
end

function overapproximate(X::Intersection{N, <:Hyperplane, <:AbstractZonotope},
                         dirs::Type{<:AbstractDirections}) where {N}
    return overapproximate(X.Y ∩ X.X, dirs(dim(X)))
end

# ===========================================================
# Functionality that requires IntervalConstraintProgramming
# ===========================================================

# function to be loaded by Requires
function load_paving_overapproximation()

return quote

using .IntervalConstraintProgramming: Paving


"""
    overapproximate(p::Paving{L, N}, dirs::AbstractDirections{N, VN}) where {L, N, VN}

Overapproximation of a Paving-type set representation using a polyhedron in constraint representation.

### Input

- `p`    -- paving
- `dirs` -- template directions

### Output

An overapproximation of a paving using a polyhedron in constraint representation (`HPolyhedron`) with constraints in direction `dirs`.

### Algorithm

This function takes the union of the elements in the boundary of p, first
converted into hyperrectangles, and then calculates the support function of the
set along each  direction in dirs, to compute the `HPolyhedron` constraints.

### Requires IntervalConstraintProgramming
"""
function overapproximate(p::Paving{L, N}, dirs::AbstractDirections{N, VN}) where {L, N, VN}
    # enclose outer approximation
    Uouter = UnionSetArray(convert.(Hyperrectangle, p.boundary))
    constraints = [HalfSpace(d, ρ(d, Uouter)) for d in dirs]
    return HPolyhedron(constraints)
end

# alias with HPolyhedron type as second argument
function overapproximate(p::Paving{L, N}, ::Type{<:HPolyhedron}, dirs::AbstractDirections{N, VN}) where {L, N, VN}
    return overapproximate(p, dirs)
end

end # quote
end # load_paving_overapproximation
