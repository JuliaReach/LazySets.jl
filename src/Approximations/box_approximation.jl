# ==================
# 2D box directions
# ==================

const DIR_EAST(N) = [one(N), zero(N)]
const DIR_NORTH(N) = [zero(N), one(N)]
const DIR_WEST(N) = [-one(N), zero(N)]
const DIR_SOUTH(N) = [zero(N), -one(N)]

dir_east(N, ::AbstractVector) = DIR_EAST(N)
dir_north(N, ::AbstractVector) = DIR_NORTH(N)
dir_west(N, ::AbstractVector) = DIR_WEST(N)
dir_south(N, ::AbstractVector) = DIR_SOUTH(N)

function load_staticarrays_directions()
    return quote
        const DIR_EAST_STATIC(N) = SVector{2}([one(N), zero(N)])
        const DIR_NORTH_STATIC(N) = SVector{2}([zero(N), one(N)])
        const DIR_WEST_STATIC(N) = SVector{2}([-one(N), zero(N)])
        const DIR_SOUTH_STATIC(N) = SVector{2}([zero(N), -one(N)])

        dir_east(N, ::SVector) = DIR_EAST_STATIC(N)
        dir_north(N, ::SVector) = DIR_NORTH_STATIC(N)
        dir_west(N, ::SVector) = DIR_WEST_STATIC(N)
        dir_south(N, ::SVector) = DIR_SOUTH_STATIC(N)
    end
end  # quote / load_staticarrays_directions

"""
    box_approximation(S::LazySet)

Overapproximate a set by a tight hyperrectangle.

### Input

- `S` -- set

### Output

A tight hyperrectangle.

### Notes

The convenience aliases `interval_hull` and `□` are also available. `□` can be
typed by `\\square<tab>`.

### Algorithm

The center and radius of the hyperrectangle are obtained by averaging the low
and high coordinates of `S` computed with the `extrema` function.
"""
function box_approximation(S::LazySet)
    return _box_approximation_extrema(S)
end

function _box_approximation_extrema(S::LazySet{N}) where {N}
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        lo, hi = extrema(S, i)
        if isinf(lo) || isinf(hi)
            throw(ArgumentError("the box approximation of an unbounded set is undefined"))
        end
        ri = (hi - lo) / 2
        if ri < zero(N)
            if !_geq(ri, zero(N))
                # contradicting bounds => set is empty
                return EmptySet{N}(n)
            else
                # assume this is a precision error with flat sets
                ri = zero(N)
            end
        end
        c[i] = (hi + lo) / 2
        r[i] = ri
    end
    return Hyperrectangle(c, r)
end

function box_approximation(cap::Intersection{N,T1,T2}) where {N,T1,T2}
    if ispolyhedral(cap.X) && ispolyhedral(cap.Y)
        if isboundedtype(T1) || isboundedtype(T2)
            S = convert(HPolytope, cap)
        else
            S = convert(HPolyhedron, cap)
        end
    else
        S = cap
    end
    return _box_approximation_extrema(S)
end

interval_hull = box_approximation

□(X::LazySet) = box_approximation(X)

# AbstractHyperrectangle specialization
function box_approximation(S::AbstractHyperrectangle)
    return Hyperrectangle(center(S), radius_hyperrectangle(S))
end

# empty set specialization
box_approximation(∅::EmptySet) = ∅

"""
    box_approximation(S::CartesianProductArray{N, <:AbstractHyperrectangle}) where {N}

Overapproximate the Cartesian product of a finite number of hyperrectangular
sets by a tight hyperrectangle.

### Input

- `S`-- Cartesian product of a finite number of hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

This method falls back to the `convert` method. Since the sets wrapped by the
Cartesian product array are hyperrectangles, this can be done without
overapproximation.
"""
function box_approximation(S::CartesianProductArray{N,<:AbstractHyperrectangle}) where {N}
    return convert(Hyperrectangle, S)
end

"""
    box_approximation(S::CartesianProduct{N, <:AbstractHyperrectangle, <:AbstractHyperrectangle}) where {N}

Overapproximate the Cartesian product of two hyperrectangular sets by a tight
hyperrectangle.

### Input

- `S`-- Cartesian product of two hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

This method falls back to the `convert` method. Since the sets wrapped by the
Cartesian product array are hyperrectangles, this can be done without
overapproximation.
"""
function box_approximation(S::CartesianProduct{N,<:AbstractHyperrectangle,<:AbstractHyperrectangle}) where {N}
    return convert(Hyperrectangle, S)
end

"""
    box_approximation(lm::LinearMap{N, <:AbstractHyperrectangle}) where {N}

Return a tight overapproximation of the linear map of a hyperrectangular set
using a hyperrectangle.

### Input

- `lm`-- linear map of a hyperrectangular set

### Output

A hyperrectangle.

### Algorithm

If `c` and `r` denote the center and vector radius of a hyperrectangle `H`,
a tight hyperrectangular overapproximation of `M * H` is obtained by
transforming `c ↦ M*c` and `r ↦ abs.(M) * r`, where `abs.(⋅)` denotes the
element-wise absolute-value operator.
"""
function box_approximation(lm::LinearMap{N,<:AbstractHyperrectangle}) where {N}
    M, X = lm.M, lm.X
    c = M * center(X)
    r = abs.(M) * radius_hyperrectangle(X)
    return Hyperrectangle(c, r)
end

"""
    box_approximation(R::Rectification{N}) where {N}

Overapproximate the rectification of a set by a tight hyperrectangle.

### Input

- `R`-- rectification of a set

### Output

A hyperrectangle.

### Algorithm

Box approximation and rectification distribute. We first check whether the
wrapped set is empty. If so, we return the empty set. Otherwise, we compute the
box approximation of the wrapped set, rectify the resulting box (which is
simple), and finally convert the resulting set to a box.
"""
function box_approximation(R::Rectification{N}) where {N}
    if isempty(R.X)
        return EmptySet{N}(dim(R))
    end
    return convert(Hyperrectangle, Rectification(box_approximation(R.X)))
end

"""
    box_approximation(Z::AbstractZonotope)

Return a tight overapproximation of a zonotope with an axis-aligned box.

### Input

- `Z` -- zonotope

### Output

A hyperrectangle.

### Algorithm

This function implements the method in [AlthoffSB10; Section 5.1.2](@citet). A zonotope
``Z = ⟨c, G⟩`` can be tightly overapproximated by an axis-aligned hyperrectangle
such that its center is ``c`` and the radius along dimension ``i`` is the
column-sum of the absolute values of the ``i``-th row of ``G`` for
``i = 1,…, p``, where ``p`` is the number of generators of ``Z``.
"""
function box_approximation(Z::AbstractZonotope)
    r = _box_radius(Z)
    return Hyperrectangle(center(Z), r)
end

"""
    box_approximation(am::AbstractAffineMap{N, <:AbstractHyperrectangle}) where {N}

Overapproximate the affine map of a hyperrectangular set by a tight
hyperrectangle.

### Input

- `am` -- affine map of a hyperrectangular set

### Output

A hyperrectangle.

### Algorithm

If `c` and `r` denote the center and vector radius of a hyperrectangle `H`
and `v` is the translation vector, a tight hyperrectangular overapproximation of
`M * H + v` is obtained by transforming `c ↦ M*c+v` and `r ↦ abs.(M) * r`, where
`abs.(⋅)` denotes the element-wise absolute-value operator.
"""
function box_approximation(am::AbstractAffineMap{N,<:AbstractHyperrectangle}) where {N}
    M, X, v = matrix(am), set(am), vector(am)
    c = M * center(X) + v
    r = abs.(M) * radius_hyperrectangle(X)
    return Hyperrectangle(c, r)
end

function box_approximation(P::Union{VPolytope,VPolygon})
    n = dim(P)
    vlist = vertices_list(P)
    @assert !isempty(vlist) "cannot overapproximate an empty polytope"

    @inbounds v1 = vlist[1]
    c = similar(v1)
    r = similar(v1)

    @inbounds for i in 1:n
        low_i = v1[i]
        high_i = v1[i]
        for v in vlist
            if v[i] > high_i
                high_i = v[i]
            elseif v[i] < low_i
                low_i = v[i]
            end
        end
        ri = (high_i - low_i) / 2
        c[i] = low_i + ri
        r[i] = ri
    end
    return Hyperrectangle(c, r)
end

# centrally-symmetric sets only require n support-function evaluations
function box_approximation(X::ACS)
    N = eltype(X)
    n = dim(X)
    c = center(X)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        r[i] = high(X, i) - c[i]
    end
    return Hyperrectangle(c, r)
end

# balls result in a hypercube with the same radius
function box_approximation(B::Union{Ball1,Ball2,BallInf,Ballp})
    return Hyperrectangle(center(B), fill(B.radius, dim(B)))
end

"""
    box_approximation(ch::ConvexHull; [algorithm]::String="box")

Overapproximate a convex hull with a tight hyperrectangle.

### Input

- `ch`        -- convex hull
- `algorithm` -- (optional; default: `"box"`) algorithm choice

### Output

A hyperrectangle.

### Algorithm

Let `X` and `Y` be the two sets of `ch`.
We make use of the following property:

```math
□(CH(X, Y))
    = □\\left( X ∪ Y \\right)
    = □\\left( □(X) ∪ □(Y) \\right)
```

If `algorithm == "extrema"`, we compute the low and high coordinates of `X` and
`Y` via `extrema`.

If `algorithm == "box"`, we instead compute the box approximations of `X` and
`Y` via `box_approximation`.

In both cases we then take the box approximation of the result.

The `"extrema"` algorithm is more efficient if `extrema` is efficient because
it does not need to allocate the intermediate hyperrectangles.
"""
function box_approximation(ch::ConvexHull; algorithm::String="box")
    if algorithm == "extrema"
        return _box_approximation_chull_extrema(first(ch), second(ch))
    elseif algorithm == "box"
        return _box_approximation_chull_box(first(ch), second(ch))
    else
        throw(ArgumentError("unknown algorithm $algorithm"))
    end
end

function _box_approximation_chull_extrema(X::LazySet{N}, Y) where {N}
    n = dim(X)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    for i in 1:n
        l1, h1 = extrema(X, i)
        l2, h2 = extrema(Y, i)
        li = min(l1, l2)
        hi = max(h1, h2)
        ci = (hi + li) / 2
        c[i] = ci
        r[i] = hi - ci
    end
    return Hyperrectangle(c, r)
end

function _box_approximation_chull_box(X::LazySet{N}, Y) where {N}
    n = dim(X)
    H1 = box_approximation(X)
    H2 = box_approximation(Y)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        li = min(low(H1, i), low(H2, i))
        hi = max(high(H1, i), high(H2, i))
        ci = (hi + li) / 2
        c[i] = ci
        r[i] = hi - ci
    end
    return Hyperrectangle(c, r)
end

"""
    box_approximation(ms::MinkowskiSum)

Overapproximate the Minkowski sum of two sets with a tight hyperrectangle.

### Input

- `ms` -- Minkowski sum

### Output

A hyperrectangle.

### Algorithm

The box approximation distributes over the Minkowski sum:

```math
□(X ⊕ Y) = □(X) ⊕ □(Y)
```

It suffices to compute the box approximation of each summand and then take the
concrete Minkowski sum for hyperrectangles.
"""
function box_approximation(ms::MinkowskiSum)
    H1 = box_approximation(first(ms))
    H2 = box_approximation(second(ms))
    return minkowski_sum(H1, H2)
end

# function to be loaded by Requires
function load_taylormodels_box_approximation()
    return quote
        using .TaylorModels: TaylorModel1, TaylorModelN, domain, evaluate

        box_approximation(vTM::Vector{<:TaylorModel1}) = _box_approximation_vTM(vTM)

        box_approximation(vTM::Vector{<:TaylorModelN}) = _box_approximation_vTM(vTM)

        function _box_approximation_vTM(vTM)
            B = IA.IntervalBox([evaluate(vTM[i], domain(p)) for (i, p) in enumerate(vTM)]...)
            return convert(Hyperrectangle, B)
        end
    end
end  # quote / load_taylormodels_box_approximation
