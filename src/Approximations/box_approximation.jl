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

using .StaticArrays: SVector

const DIR_EAST_STATIC(N) = SVector{2}([one(N), zero(N)])
const DIR_NORTH_STATIC(N) = SVector{2}([zero(N), one(N)])
const DIR_WEST_STATIC(N) = SVector{2}([-one(N), zero(N)])
const DIR_SOUTH_STATIC(N) = SVector{2}([zero(N), -one(N)])

dir_east(N, ::SVector) = DIR_EAST_STATIC(N)
dir_north(N, ::SVector) = DIR_NORTH_STATIC(N)
dir_west(N, ::SVector) = DIR_WEST_STATIC(N)
dir_south(N, ::SVector) = DIR_SOUTH_STATIC(N)

end end  # quote / load_staticarrays_directions

# ================================================
# Concrete overapproximation with a hyperrectangle
# ================================================

"""
    box_approximation(S::LazySet{N}) where {N}

Overapproximate a set by a tight hyperrectangle.

### Input

- `S` -- set

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet{N}) where {N}
    c, r = box_approximation_helper(S)
    if r[1] < 0
        return EmptySet{N}(dim(S))
    end
    return Hyperrectangle(c, r)
end

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    box_approximation_helper(S::LazySet{N}) where {N}

Common code of `box_approximation` and `box_approximation_symmetric`.

### Input

- `S` -- set

### Output

A tuple containing the data that is needed to construct a tightly
overapproximating hyperrectangle.

- `c` -- center
- `r` -- radius

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.
"""
@inline function box_approximation_helper(S::LazySet{N}) where {N}
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        htop = ρ(SingleEntryVector(i, n, one(N)), S)
        hbottom = -ρ(SingleEntryVector(i, n, -one(N)), S)
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
        if r[i] < 0
            # contradicting bounds => set is empty
            # terminate with first radius entry being negative
            r[1] = r[i]
            break
        end
    end
    return c, r
end

# ===============
# Specializations
# ===============

# hyperrectangle specialization
function box_approximation(S::AbstractHyperrectangle)
    return Hyperrectangle(center(S), radius_hyperrectangle(S))
end

# empty set specialization
box_approximation(∅::EmptySet) = ∅

"""
    box_approximation(S::CartesianProductArray{N, <:AbstractHyperrectangle}) where {N}

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
function box_approximation(S::CartesianProductArray{N, <:AbstractHyperrectangle}) where {N}
    return convert(Hyperrectangle, S)
end

"""
    box_approximation(S::CartesianProduct{N, <:AbstractHyperrectangle, <:AbstractHyperrectangle}) where {N}

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
function box_approximation(S::CartesianProduct{N, <:AbstractHyperrectangle, <:AbstractHyperrectangle}) where {N}
    return convert(Hyperrectangle, S)
end

"""
    box_approximation(lm::LinearMap{N, <:AbstractHyperrectangle}) where {N}

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
`c ↦ M*c` and `r ↦ abs.(M) * r`, where `abs.(⋅)` denotes the element-wise absolute
value operator.
"""
function box_approximation(lm::LinearMap{N, <:AbstractHyperrectangle}) where {N}
    M, X = lm.M, lm.X
    center_MX = M * center(X)
    radius_MX = abs.(M) * radius_hyperrectangle(X)
    return Hyperrectangle(center_MX, radius_MX)
end

"""
    box_approximation(r::Rectification{N}) where {N}

Overapproximate the rectification of a convex set by a tight hyperrectangle.

### Input

- `r`              -- rectification of a convex set
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

Box approximation and rectification distribute.
Hence we first check whether the wrapped set is empty.
If so, we return the empty set.
Otherwise, we compute the box approximation of the wrapped set, rectify the
resulting box (which is simple), and finally convert the resulting set to a box.
"""
function box_approximation(r::Rectification{N}) where {N}
    if isempty(r.X)
        return EmptySet{N}(dim(r))
    end
    return convert(Hyperrectangle, Rectification(box_approximation(r.X)))
end

"""
    box_approximation(Z::AbstractZonotope)

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
function box_approximation(Z::AbstractZonotope)
    r = sum(abs.(genmat(Z)), dims=2)[:]
    return Hyperrectangle(center(Z), r)
end

"""
    box_approximation(am::AbstractAffineMap{N, <:AbstractHyperrectangle}) where {N}

Overapproximate the affine map of a hyperrectangular set
using a hyperrectangle.

### Input

- `am`             -- affine map of a hyperrectangular set
- `Hyperrectangle` -- type for dispatch

### Output

A hyperrectangle.

### Algorithm

If `c` and `r` denote the center and vector radius of a hyperrectangle `H`
and `v` the translation vector, a tight hyperrectangular overapproximation of
`M * H + v` is obtained by transforming `c ↦ M*c+v` and `r ↦ abs.(M) * r`, where
`abs.(⋅)` denotes the element-wise absolute value operator.
"""
function box_approximation(am::AbstractAffineMap{N, <:AbstractHyperrectangle}) where {N}
    M, X, v = matrix(am), set(am), vector(am)
    center_MXv = M * center(X) + v
    radius_MX = abs.(M) * radius_hyperrectangle(X)
    return Hyperrectangle(center_MXv, radius_MX)
end

function box_approximation(P::Union{VPolytope, VPolygon})
    n = dim(P)
    vlist = vertices_list(P)
    @assert !isempty(vlist) "cannot overapproximate an empty polytope"

    @inbounds v1 = vlist[1]
    center = similar(v1)
    radius = similar(v1)

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
        radius_i = (high_i - low_i) / 2
        center[i] = low_i + radius_i
        radius[i] = radius_i
    end
    return Hyperrectangle(center, radius)
end

# centrally symmetric sets only require n support-function evaluations
function box_approximation(X::AbstractCentrallySymmetric{N}) where {N}
    n = dim(X)
    c = center(X)
    r = Vector{N}(undef, n)
    @inbounds for i in 1:n
        r[i] = ρ(SingleEntryVector(i, n, one(N)), X) - c[i]
    end
    return Hyperrectangle(c, r)
end

# balls result in a hypercube with the same radius
function box_approximation(B::Union{Ball1, Ball2, BallInf, Ballp})
    return Hyperrectangle(center(B), fill(B.radius, dim(B)))
end
