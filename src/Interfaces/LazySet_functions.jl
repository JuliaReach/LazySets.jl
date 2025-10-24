export neutral,
       absorbing,
       tosimplehrep,
       singleton_list,
       chebyshev_center_radius,
       flatten,
       triangulate,
       triangulate_faces

"""
    triangulate(X::LazySet; [algorithm]::String="delaunay", [kwargs]...)

Compute a triangulation of the given polytopic set.

### Input

- `X`         -- polytopic set
- `algorithm` -- (optional; default: `delaunay`) string to choose the
                 type of triangulation
- `kwargs`    -- further keyword arguments passed to the algorithm

### Output

A union of polytopes in vertex representation.

### Algorithm

The algorithm is selected with the argument `algorithm`.

- `"delaunay"`: This algorithm computes a Delaunay triangulation.

This algorithm can receive another optional argument `compute_triangles_3d`,
a Boolean flag that defaults to `false`. It is used to compute the 2D
triangulation of a 3D set if `true`.

The implementation requires the package
[MiniQhull.jl](https://github.com/gridap/MiniQhull.jl), which uses the library
[Qhull](http://www.qhull.org/).

The algorithm works in arbitrary dimension and requires that the list of
vertices of `X` can be obtained.
"""
function triangulate(X::LazySet; algorithm::String="delaunay", kwargs...)
    if algorithm == "delaunay"
        require(@__MODULE__, :MiniQhull; fun_name="triangulate")
        return _triangulate_delaunay(X; kwargs...)
    else
        throw(ArgumentError("unknown algorithm $algorithm"))
    end
end

@validate function triangulate_faces(X)
    require(@__MODULE__, :Polyhedra; fun_name="triangulate_faces")
    throw(ArgumentError("`triangulate_faces` not implemented for $(typeof(X))"))
end

"""
# Extended help

    isconvextype(::Type{<:LazySet})

### Algorithm

The default implementation returns `false`. All set types that can determine
convexity should override this behavior.

### Examples

A ball in the infinity norm is always convex:

```jldoctest convex_types
julia> isconvextype(BallInf)
true
```

The union (`UnionSet`) of two sets may in general be either convex or not. Since
convexity cannot be decided by just using type information, `isconvextype`
returns `false`.

```jldoctest convex_types
julia> isconvextype(UnionSet)
false
```

However, the type parameters of set operations allow to decide convexity in some
cases by falling back to the convexity information of the argument types.
Consider the lazy intersection. The intersection of two convex sets is always
convex:

```jldoctest convex_types
julia> isconvextype(Intersection{Float64, BallInf{Float64}, BallInf{Float64}})
true
```
"""
isconvextype(::Type{<:LazySet}) = false

# Polyhedra backend (fallback method)
function default_polyhedra_backend(P::LazySet{N}) where {N}
    require(@__MODULE__, :Polyhedra; fun_name="default_polyhedra_backend")
    if dim(P) == 1
        return default_polyhedra_backend_1d(N)
    else
        return default_polyhedra_backend_nd(N)
    end
end

# Note: this method cannot be documented due to a bug in Julia
@validate function low(X::LazySet, i::Int)
    return _low(X, i)
end

function _low(X::LazySet{N}, i::Int) where {N}
    n = dim(X)
    d = SingleEntryVector(i, n, -one(N))
    return -ρ(d, X)
end

function _low_vlist(X::LazySet, i::Int)
    return _low_vlist(vertices_list(X), i)
end

function _low_vlist(vlist::AbstractVector{<:AbstractVector{N}}, i::Int) where {N}
    l = N(Inf)
    @inbounds for v in vlist
        l = min(l, v[i])
    end
    return l
end

"""
# Extended help

    low(X::LazySet)

### Algorithm

The default implementation applies `low` in each dimension.
"""
function low(X::LazySet)
    return _low(X)
end

function _low(X::LazySet)
    n = dim(X)
    return [low(X, i) for i in 1:n]
end

function _low_vlist(X::LazySet)
    n = dim(X)
    vlist = vertices_list(X)
    return [_low_vlist(vlist, i) for i in 1:n]
end

# Note: this method cannot be documented due to a bug in Julia
@validate function high(X::LazySet, i::Int)
    return _high(X, i)
end

function _high(X::LazySet{N}, i::Int) where {N}
    n = dim(X)
    d = SingleEntryVector(i, n, one(N))
    return ρ(d, X)
end

function _high_vlist(X::LazySet, i::Int)
    return _high_vlist(vertices_list(X), i)
end

function _high_vlist(vlist::AbstractVector{<:AbstractVector{N}}, i::Int) where {N}
    h = N(-Inf)
    @inbounds for v in vlist
        h = max(h, v[i])
    end
    return h
end

"""
# Extended help

    high(X::LazySet)

### Algorithm

The default implementation applies `high` in each dimension.
"""
function high(X::LazySet)
    return _high(X)
end

function _high(X::LazySet)
    n = dim(X)
    return [high(X, i) for i in 1:n]
end

function _high_vlist(X::LazySet)
    n = dim(X)
    vlist = vertices_list(X)
    return [_high_vlist(vlist, i) for i in 1:n]
end

"""
# Extended help

    extrema(X::LazySet, i::Int)

### Algorithm

The default implementation computes the extrema via `low` and `high`.
"""
@validate function extrema(X::LazySet, i::Int)
    return _extrema_lowhigh(X, i)
end

function _extrema_lowhigh(X::LazySet, i::Int)
    l = low(X, i)
    h = high(X, i)
    return (l, h)
end

function _extrema_vlist(X::LazySet, i::Int)
    return _extrema_vlist(vertices_list(X), i)
end

function _extrema_vlist(vlist::AbstractVector{<:AbstractVector{N}}, i::Int) where {N}
    if isempty(vlist)
        return (N(Inf), N(-Inf))
    end
    l = h = @inbounds vlist[1][i]
    @inbounds for v in vlist
        vi = v[i]
        if vi > h
            h = vi
        elseif vi < l
            l = vi
        end
    end
    return (l, h)
end

"""
# Extended help

    extrema(X::LazySet)

### Algorithm

The default implementation computes the extrema via `low` and `high`.
"""
function extrema(X::LazySet)
    return _extrema_lowhigh(X)
end

function _extrema_lowhigh(X::LazySet)
    l = low(X)
    h = high(X)
    return (l, h)
end

function _extrema_vlist(X::LazySet{N}) where {N}
    n = dim(X)
    if n <= 0
        return (N[], N[])
    end
    vlist = vertices_list(X)
    if isempty(vlist)
        return (fill(N(Inf), n), fill(N(-Inf), n))
    end
    l = copy(@inbounds vlist[1])
    h = copy(l)
    @inbounds for v in vlist
        for i in 1:n
            vi = v[i]
            if vi > h[i]
                h[i] = vi
            elseif vi < l[i]
                l[i] = vi
            end
        end
    end
    return (l, h)
end

"""
    plot_vlist(X::S, ε::Real) where {S<:LazySet}

Return a list of vertices used for plotting a two-dimensional set.

### Input

- `X` -- two-dimensional set
- `ε` -- precision parameter

### Output

A list of vertices of a polygon `P`.
For convex `X`, `P` usually satisfies that the Hausdorff distance to `X` is less
than `ε`.
"""
function plot_vlist(X::S, ε::Real) where {S<:LazySet}
    @assert isconvextype(S) "can only plot convex sets"

    P = overapproximate(X, ε)
    return convex_hull(vertices_list(P))
end

"""
# Extended help

    convex_hull(X::LazySet; kwargs...)

### Output

The set `X` itself if its type indicates that it is convex, or a new set with
the list of the vertices describing the convex hull.

### Algorithm

For non-convex sets, this method relies on the `vertices_list` method.
"""
function convex_hull(X::LazySet; kwargs...)
    if isconvextype(typeof(X))
        return X
    end
    return _convex_hull_polytopes(X; kwargs...)
end

function _convex_hull_polytopes(X; kwargs...)
    vlist = convex_hull(vertices_list(X); kwargs...)
    return _convex_hull_set(vlist; n=dim(X))
end

"""
# Extended help

    eltype(::Type{<:LazySet{N}}) where {N}

### Algorithm

The default implementation assumes that the first type parameter is the numeric
type.
"""
eltype(::Type{<:LazySet{N}}) where {N} = N

"""
# Extended help

    eltype(::LazySet{N}) where {N}

### Algorithm

The default implementation assumes that the first type parameter is the numeric
type.
"""
eltype(::LazySet{N}) where {N} = N

"""
# Extended help

    ρ(d::AbstractVector, X::LazySet)

### Algorithm

The default implementation computes a support vector via `σ`.
"""
@validate function ρ(d::AbstractVector, X::LazySet)
    return dot(d, σ(d, X))
end

"""
# Extended help

    isboundedtype(::Type{<:LazySet})

### Notes

Note that some sets may still represent an unbounded set even though their type
actually does not (example: [`HPolytope`](@ref), because the construction with
non-bounding linear constraints is allowed).

### Algorithm

The default implementation returns `false`. All set types that can determine
boundedness should override this behavior.
"""
function isboundedtype(::Type{<:LazySet})
    return false
end

"""
# Extended help

    isbounded(X::LazySet; [algorithm]="support_function")

### Input

- `algorithm` -- (optional, default: `"support_function"`) algorithm choice,
                 possible options are `"support_function"` and `"stiemke"`

### Algorithm

See the documentation of `_isbounded_unit_dimensions` or `_isbounded_stiemke`
for details.
"""
function isbounded(X::LazySet; algorithm="support_function")
    if algorithm == "support_function"
        return _isbounded_unit_dimensions(X)
    elseif algorithm == "stiemke"
        return _isbounded_stiemke(constraints_list(X))
    else
        throw(ArgumentError("unknown algorithm $algorithm"))
    end
end

"""
    _isbounded_unit_dimensions(X::LazySet)

Check whether a set is bounded in each unit dimension.

### Input

- `X` -- set

### Output

`true` iff the set is bounded in each unit dimension.

### Algorithm

This function asks for upper and lower bounds in each ambient dimension.
"""
function _isbounded_unit_dimensions(X::LazySet)
    @inbounds for i in 1:dim(X)
        if isinf(low(X, i)) || isinf(high(X, i))
            return false
        end
    end
    return true
end

"""
# Extended help

    ispolyhedraltype(T::Type{<:LazySet})

### Algorithm

The default implementation returns `false`.
"""
function ispolyhedraltype(::Type{<:LazySet})
    return false
end

"""
# Extended help

    ispolyhedral(X::LazySet)

### Algorithm

The default implementation checks `ispolyhedraltype(typeof(X))`.
"""
function ispolyhedral(X::LazySet)
    return ispolyhedraltype(typeof(X))
end

"""
# Extended help

    ispolytopic(X::LazySet)

### Algorithm

The default implementation checks `ispolyhedral(X)` and `isbounded(X)`. This is
typically enough, but note that these functions may give a conservative result.
"""
function ispolytopic(X::LazySet)
    return ispolyhedral(X) && isbounded(X)
end

"""
# Extended help

    norm(X::LazySet, [p]::Real=Inf)

### Algorithm

The default implementation handles the case `p == Inf` using `extrema`.
Otherwise it checks whether `X` is polytopic, in which case it iterates over all
vertices.
"""
@validate function norm(X::LazySet, p::Real=Inf)
    return _norm_default(X, p)
end

function _norm_default(X::LazySet, p::Real)
    if p == Inf
        l, h = extrema(X)
        return max(maximum(abs, l), maximum(abs, h))
    elseif ispolytopic(X)
        return maximum(norm(v, p) for v in vertices_list(X))
    else
        error("the norm for this value of p=$p is not implemented")
    end
end

"""
# Extended help

    radius(X::LazySet, [p]::Real=Inf)

### Algorithm

The default implementation handles the case `p == Inf` using
`ballinf_approximation`.
"""
@validate function radius(X::LazySet, p::Real=Inf)
    if p == Inf
        return radius(Approximations.ballinf_approximation(X), p)
    else
        error("the radius for this value of p=$p is not implemented")
    end
end

"""
# Extended help

    diameter(X::LazySet, [p]::Real=Inf)

### Algorithm

The default implementation applies the function `radius` and doubles the result.
"""
@validate function diameter(X::LazySet, p::Real=Inf)
    return radius(X, p) * 2
end

"""
# Extended help

    translate(X::LazySet, v::AbstractVector)

### Algorithm

The default implementation calls `translate!` on a copy of `X`.
"""
@validate function translate(X::LazySet, v::AbstractVector)
    Y = copy(X)
    translate!(Y, v)
    return Y
end

"""
# Extended help

    affine_map(M, X::LazySet, v::AbstractVector; [kwargs]...)

### Algorithm

The default implementation applies the functions `linear_map` and `translate`.
"""
@validate function affine_map(M, X::LazySet, v::AbstractVector; kwargs...)
    return translate(linear_map(M, X; kwargs...), v)
end

"""
# Extended help

    exponential_map(M::AbstractMatrix, X::LazySet)

### Algorithm

The default implementation applies the functions `exp` and `linear_map`.
"""
@validate function exponential_map(M::AbstractMatrix, X::LazySet)
    return linear_map(exp(M), X)
end

"""
# Extended help

    an_element(X::LazySet)

### Algorithm

The default implementation computes a support vector along direction
``[1, 0, …, 0]``. This may fail for unbounded sets.
"""
@validate function an_element(X::LazySet)
    return _an_element_lazyset(X)
end

function _an_element_lazyset(X::LazySet)
    N = eltype(X)
    e₁ = SingleEntryVector(1, dim(X), one(N))
    v = σ(e₁, X)
    if any(isinf, v)
        throw(ArgumentError("this implementation assumes a bounded set"))
    end
    return v
end

# hook into random API
function rand(rng::AbstractRNG, ::SamplerType{T}) where {T<:LazySet}
    return rand(T; rng=rng)
end

# This function computes a `copy` of each field in `X`. See the documentation of
# `?copy` for further details.
function copy(X::T) where {T<:LazySet}
    args = [copy(getfield(X, f)) for f in fieldnames(T)]
    BT = basetype(X)
    return BT(args...)
end

"""
    tosimplehrep(S::LazySet)

Return the simple constraint representation ``Ax ≤ b`` of a polyhedral set from
its list of linear constraints.

### Input

- `S` -- polyhedral set

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.

### Algorithm

This fallback implementation relies on `constraints_list(S)`.
"""
tosimplehrep(S::LazySet) = tosimplehrep(constraints_list(S))

"""
# Extended help

    reflect(P::LazySet)

### Output

The set `-P`, which is either of type `HPolytope` if `P` is a polytope (i.e.,
bounded) or of type `HPolyhedron` otherwise.

### Algorithm

This function requires that the list of constraints of the set `P` is
available, i.e., that it can be written as
``P = \\{z ∈ ℝⁿ: ⋂ sᵢᵀz ≤ rᵢ, i = 1, ..., N\\}.``

This function can be used to implement the alternative definition of the
Minkowski difference
```math
    A ⊖ B = \\{a − b \\mid a ∈ A, b ∈ B\\} = A ⊕ (-B)
```
by calling `minkowski_sum(A, reflect(B))`.
"""
function reflect(P::LazySet)
    if !ispolyhedral(P)
        error("this implementation requires a polyhedral set; try " *
              "overapproximating with an `HPolyhedron` first")
    end

    F, g = tosimplehrep(P)
    T = isbounded(P) ? HPolytope : HPolyhedron
    return T(-F, g)
end

function is_interior_point(v::AbstractVector{<:Real}, X::LazySet; kwargs...)
    N = promote_type(eltype(v), eltype(X))
    if N != eltype(X)
        throw(ArgumentError("the set eltype must be more general"))
    end
    if N != eltype(v)
        v = convert(Vector{N}, v)
    end
    ε = get(kwargs, :ε, _rtol(N))
    p = get(kwargs, :p, N(Inf))
    return is_interior_point(v, X; p=p, ε=ε)
end

"""
# Extended help

    is_interior_point(v::AbstractVector{N}, X::LazySet{N}; [p]=N(Inf), [ε]=_rtol(N)) where {N<:Real}

### Algorithm

The default implementation determines `v ∈ interior(X)` with error tolerance
`ε` by checking whether a `Ballp` of norm `p` with center `v` and radius `ε` is
contained in `X`.
"""
@validate function is_interior_point(v::AbstractVector{N}, X::LazySet{N}; p=N(Inf),
                                     ε=_rtol(N)) where {N<:Real}
    return Ballp(p, v, ε) ⊆ X
end

"""
    plot_recipe(X::LazySet, [ε])

Convert a compact convex set to a pair `(x, y)` of points for plotting.

### Input

- `X` -- compact convex set
- `ε` -- approximation-error bound

### Output

A pair `(x, y)` of points that can be plotted.

### Notes

We do not support three-dimensional or higher-dimensional sets at the moment.

### Algorithm

One-dimensional sets are converted to an `Interval`.

For two-dimensional sets, we first compute a polygonal overapproximation.
The second argument, `ε`, corresponds to the error in Hausdorff distance between
the overapproximating set and `X`.
On the other hand, if you only want to produce a fast box-overapproximation of
`X`, pass `ε=Inf`.

Finally, we use the plot recipe for the constructed set (interval or polygon).
"""
function plot_recipe(X::LazySet, ε)
    @assert dim(X) <= 3 "cannot plot a $(dim(X))-dimensional $(typeof(X))"
    @assert isboundedtype(typeof(X)) || isbounded(X) "cannot plot an " *
                                                     "unbounded $(typeof(X))"
    @assert isconvextype(typeof(X)) "can only plot convex sets"

    if dim(X) == 1
        Y = convert(Interval, X)
    elseif dim(X) == 2
        Y = overapproximate(X, ε)
    else
        return _plot_recipe_3d_polytope(X)
    end
    return plot_recipe(Y, ε)
end

function _plot_recipe_3d_polytope(P::LazySet, N=eltype(P))
    require(@__MODULE__, :MiniQhull; fun_name="_plot_recipe_3d_polytope")
    @assert ispolytopic(P) "3D plotting is only available for polytopes"

    vlist, C = triangulate_vlist_connectivity(P; compute_triangles_3d=true)

    m = length(vlist)
    if m == 0
        @warn "received a polyhedron with no vertices during plotting"
        return plot_recipe(EmptySet{N}(2), zero(N))
    end

    x = Vector{N}(undef, m)
    y = Vector{N}(undef, m)
    z = Vector{N}(undef, m)
    @inbounds for (i, vi) in enumerate(vlist)
        x[i] = vi[1]
        y[i] = vi[2]
        z[i] = vi[3]
    end

    l = size(C, 2)
    i = Vector{Int}(undef, l)
    j = Vector{Int}(undef, l)
    k = Vector{Int}(undef, l)
    @inbounds for idx in 1:l
        # normalization: -1 for zero indexing; convert to Int on 64-bit systems
        i[idx] = Int(C[1, idx] - 1)
        j[idx] = Int(C[2, idx] - 1)
        k[idx] = Int(C[3, idx] - 1)
    end

    return x, y, z, i, j, k
end

"""
# Extended help

    isoperation(X::LazySet)

### Algorithm

The default implementation checks whether the set type of the input is an
operation type using [`isoperationtype(::Type{<:LazySet})`](@ref).

### Examples

```jldoctest
julia> B = BallInf([0.0, 0.0], 1.0);

julia> isoperation(B)
false

julia> isoperation(B ⊕ B)
true
```
"""
function isoperation(X::LazySet)
    return isoperationtype(typeof(X))
end

# common error
function isoperation(::Type{<:LazySet})
    throw(ArgumentError("`isoperation` cannot be applied to a set type; use " *
                        "`isoperationtype` instead"))
end

# common error
function isoperationtype(::LazySet)
    throw(ArgumentError("`isoperationtype` cannot be applied to a set " *
                        "instance; use `isoperation` instead"))
end

"""
# Extended help

    area(X::LazySet)

### Notes

This algorithm is applicable to any two-dimensional polytopic set `X` whose list
of vertices can be computed via `vertices_list`.

### Algorithm

Let `m` be the number of vertices of `X`. We consider the following instances:

- `m = 0, 1, 2`: the output is zero.
- `m = 3`: the triangle case is solved using the Shoelace formula with 3 points.
- `m = 4`: the quadrilateral case is solved by the factored version of the
           Shoelace formula with 4 points.

Otherwise, the general Shoelace formula is used; for details see the
[Wikipedia page](https://en.wikipedia.org/wiki/Shoelace_formula).
"""
@validate function area(X::LazySet)
    @assert dim(X) == 2 "this implementation only applies to two-dimensional sets, " *
                        "but the given set is $(dim(X))-dimensional"
    @assert ispolytopic(X) "this method requires a polytope"

    vlist = vertices_list(X)
    return _area_vlist_2D(vlist)
end

# Notes:
# - dimension is expected to be 2D
# - implementation requires sorting of vertices
# - convex hull is applied in-place
function _area_vlist_2D(vlist; apply_convex_hull::Bool=true)
    if apply_convex_hull
        convex_hull!(vlist)
    end

    m = length(vlist)

    if m <= 2
        N = eltype(eltype(vlist))
        return zero(N)
    end

    if m == 3 # triangle
        res = _area_triangle(vlist)

    elseif m == 4 # quadrilateral
        res = _area_quadrilateral(vlist)

    else # general case
        res = _area_polygon(vlist)
    end

    return res
end

function _area_triangle(v::Vector{VN}) where {N,VN<:AbstractVector{N}}
    A = v[1]
    B = v[2]
    C = v[3]
    res = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2])
    return abs(res / 2)
end

function _area_quadrilateral(v::Vector{VN}) where {N,VN<:AbstractVector{N}}
    A = v[1]
    B = v[2]
    C = v[3]
    D = v[4]
    res = A[1] * (B[2] - D[2]) + B[1] * (C[2] - A[2]) + C[1] * (D[2] - B[2]) +
          D[1] * (A[2] - C[2])
    return abs(res / 2)
end

function _area_polygon(v::Vector{VN}) where {N,VN<:AbstractVector{N}}
    m = length(v)
    @inbounds res = v[m][1] * v[1][2] - v[1][1] * v[m][2]
    for i in 1:(m - 1)
        @inbounds res += v[i][1] * v[i + 1][2] - v[i + 1][1] * v[i][2]
    end
    return abs(res / 2)
end

"""
    singleton_list(P::LazySet)

Return the vertices of a polytopic set as a list of singletons.

### Input

- `P` -- polytopic set

### Output

A list of the vertices of `P` as `Singleton`s.

### Notes

This function relies on `vertices_list`, which raises an error if the set is
not polytopic (e.g., unbounded).
"""
function singleton_list(P::LazySet)
    return [Singleton(x) for x in vertices_list(P)]
end

"""
# Extended help

    concretize(X::LazySet)

### Algorithm

The default implementation returns `X`. All relevant lazy set types should
override this behavior, typically by recursively calling `concretize` on the
set arguments.
"""
function concretize(X::LazySet)
    if isoperationtype(typeof(X))
        throw(ArgumentError("concretizing $(typeof(X)) is not implemented"))
    end
    return X
end

"""
# Extended help

    constraints(X::LazySet)

### Algorithm

The default implementation computes all constraints via `constraints_list`.
"""
function constraints(X::LazySet)
    return constraints_list(X)
end

"""
# Extended help

    vertices(X::LazySet)

### Algorithm

The default implementation computes all vertices via `vertices_list`.
"""
function vertices(X::LazySet)
    return vertices_list(X)
end

function load_MiniQhull_triangulate()
    return quote
        function _triangulate_delaunay(X::LazySet; compute_triangles_3d::Bool=false)
            vlist, connect_mat = _delaunay_vlist_connectivity(X;
                                                              compute_triangles_3d=compute_triangles_3d)
            nsimplices = size(connect_mat, 2)
            if compute_triangles_3d
                simplices = [VPolytope(vlist[connect_mat[1:3, j]]) for j in 1:nsimplices]
            else
                simplices = [VPolytope(vlist[connect_mat[:, j]]) for j in 1:nsimplices]
            end
            return UnionSetArray(simplices)
        end

        # compute the vertices and the connectivity matrix of the Delaunay triangulation
        #
        # if the flag `compute_triangles_3d` is set, the resulting matrix still has four
        # rows, but the last row has no meaning
        function _delaunay_vlist_connectivity(X::LazySet;
                                              compute_triangles_3d::Bool=false)
            n = dim(X)
            @assert !compute_triangles_3d || n == 3 "the `compute_triangles_3d` " *
                                                    "option requires 3D inputs"
            vlist = vertices_list(X)
            m = length(vlist)
            coordinates = vcat(vlist...)
            flags = compute_triangles_3d ? "qhull Qt" : nothing
            connectivity_matrix = MiniQhull.delaunay(n, m, coordinates, flags)
            return vlist, connectivity_matrix
        end
    end
end  # load_MiniQhull_triangulate

"""
# Extended help

    complement(X::LazySet)

### Algorithm

The default implementation assumes that `X` is polyhedral and returns a
`UnionSetArray` of `HalfSpace`s, i.e., the output is the union of the linear
constraints which are obtained by complementing each constraint of `X`. For any
pair of sets ``(X, Y)`` we have the identity ``(X ∩ Y)^C = X^C ∪ Y^C``. We can
apply this identity for each constraint that defines a polyhedral set.
"""
function complement(X::LazySet)
    if !ispolyhedral(X)
        throw(ArgumentError("this implementation requires a polyhedral set"))
    end
    return UnionSetArray(constraints_list(Complement(X)))
end

"""
    project(X::LazySet, block::AbstractVector{Int}, [::Nothing=nothing],
            [n]::Int=dim(X); [kwargs...])

Project a set to a given block by using a concrete linear map.

### Input

- `X`       -- set
- `block`   -- block structure - a vector with the dimensions of interest
- `nothing` -- (default: `nothing`) needed for dispatch
- `n`       -- (optional, default: `dim(X)`) ambient dimension of the set `X`

### Output

A set representing the projection of `X` to block `block`.

### Algorithm

We apply the function `linear_map`.
"""
@inline function project(X::LazySet, block::AbstractVector{Int},
                         ::Nothing=nothing, n::Int=dim(X); kwargs...)
    return _project_linear_map(X, block, n; kwargs...)
end

@inline function _project_linear_map(X::LazySet{N}, block::AbstractVector{Int},
                                     n::Int=dim(X); kwargs...) where {N}
    M = projection_matrix(block, n, N)
    return linear_map(M, X)
end

"""
    project(X::LazySet, block::AbstractVector{Int}, set_type::Type{<:LazySet},
            [n]::Int=dim(X); [kwargs...])

Project a set to a given block and set type, possibly involving an
overapproximation.

### Input

- `X`        -- set
- `block`    -- block structure - a vector with the dimensions of interest
- `set_type` -- target set type
- `n`        -- (optional, default: `dim(X)`) ambient dimension of the set `X`

### Output

A set of type `set_type` representing an overapproximation of the projection of
`X`.

### Algorithm

1. Project the set `X` with `M⋅X`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected set using `overapproximate` and `set_type`.
"""
@inline function project(X::LazySet, block::AbstractVector{Int},
                         set_type::Type{<:LazySet}, n::Int=dim(X);
                         kwargs...)
    lm = project(X, block, LinearMap, n)
    return overapproximate(lm, set_type)
end

"""
    project(X::LazySet, block::AbstractVector{Int},
            set_type_and_precision::Pair{<:UnionAll,<:Real}, [n]::Int=dim(X);
            [kwargs...])

Project a set to a given block and set type with a certified error bound.

### Input

- `X`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type_and_precision` -- pair `(T, ε)` of a target set type `T` and an
                              error bound `ε` for approximation
- `n`     -- (optional, default: `dim(X)`) ambient dimension of the set `X`

### Output

A set representing the epsilon-close approximation of the projection of `X`.

### Notes

Currently we only support `HPolygon` as set type, which implies that the set
must be two-dimensional.

### Algorithm

1. Project the set `X` with `M⋅X`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected set with the given error bound `ε`.
"""
@inline function project(X::LazySet, block::AbstractVector{Int},
                         set_type_and_precision::Pair{<:UnionAll,<:Real}, n::Int=dim(X);
                         kwargs...)
    set_type, ε = set_type_and_precision
    @assert length(block) == 2 && set_type == HPolygon "currently only 2D HPolygon " *
                                                       "projection is supported"

    lm = project(X, block, LinearMap, n)
    return overapproximate(lm, set_type, ε)
end

"""
    project(X::LazySet, block::AbstractVector{Int}, ε::Real, [n]::Int=dim(X);
            [kwargs...])

Project a set to a given block and set type with an error bound.

### Input

- `X`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `ε`     -- error bound for approximation
- `n`     -- (optional, default: `dim(X)`) ambient dimension of the set `X`

### Output

A set representing the epsilon-close approximation of the projection of `X`.

### Algorithm

1. Project the set `X` with `M⋅X`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected set with the given error bound `ε`.
The target set type is chosen automatically.
"""
@inline function project(X::LazySet, block::AbstractVector{Int}, ε::Real,
                         n::Int=dim(X); kwargs...)
    # currently we only support HPolygon
    if length(block) == 2
        set_type = HPolygon
    else
        throw(ArgumentError("ε-close approximation is only supported for 2D blocks"))
    end
    return project(X, block, set_type => ε, n)
end

"""
# Extended help

    rectify(X::LazySet, [concrete_intersection]::Bool=false)

### Algorithm

For each dimension in which `X` is both positive and negative, we split `X` into
these two parts and then project the negative part to zero.
"""
function rectify(X::LazySet, concrete_intersection::Bool=false)
    return to_union_of_projections(Rectification(X), concrete_intersection)
end

"""
    rationalize(::Type{T}, X::LazySet, tol::Real) where {T<:Integer}

Approximate a set represented with floating-point numbers as a set represented
with rationals of the given integer type.

### Input

- `T`   -- (optional, default: `Int`) integer type to represent the rationals
- `X`   -- set represented by floating-point components
- `tol` -- (optional, default: `eps(N)`) tolerance of the result; each rationalized
           component will differ by no more than `tol` with respect to the floating-point value

### Output

A set of the same base type as `X` where each numerical component is of type
`Rational{T}`.

### Algorithm

The default implementation rationalizes each field.
"""
function rationalize(::Type{T}, X::LazySet{<:AbstractFloat}, tol::Real) where {T<:Integer}
    m = length(fieldnames(typeof(X)))
    frat = ntuple(fi -> rationalize(T, getfield(X, fi), tol), m)
    ST = basetype(X)
    return ST(frat...)
end

# no integer type specified
rationalize(X::LazySet{<:AbstractFloat}; kwargs...) = rationalize(Int, X; kwargs...)

# `tol` as kwarg
function rationalize(::Type{T}, X::LazySet{N}; tol::Real=eps(N)) where {T<:Integer,N<:AbstractFloat}
    return rationalize(T, X, tol)
end

# vectors of sets
function rationalize(::Type{T}, X::AbstractVector{<:LazySet{<:AbstractFloat}},
                     tol::Real) where {T<:Integer}
    return rationalize.(Ref(T), X, Ref(tol))
end

"""
    chebyshev_center_radius(P::LazySet;
                            [backend]=default_polyhedra_backend(P),
                            [solver]=default_lp_solver_polyhedra(eltype(P); presolve=true))

Compute a [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of a polytopic set.

### Input

- `P`       -- polytopic set
- `backend` -- (optional; default: `default_polyhedra_backend(P)`) the backend
               for polyhedral computations
- `solver`  -- (optional; default:
               `default_lp_solver_polyhedra(N; presolve=true)`) the LP solver
               passed to `Polyhedra`

### Output

The pair `(c, r)` where `c` is a Chebyshev center of `P` and `r` is the radius
of the largest Euclidean ball with center `c` enclosed by `P`.

### Notes

The Chebyshev center is the center of a largest Euclidean ball enclosed by `P`.
In general, the center of such a ball is not unique, but the radius is.

### Algorithm

We call `Polyhedra.chebyshevcenter`.
"""
function chebyshev_center_radius(P::LazySet;
                                 backend=default_polyhedra_backend(P),
                                 solver=default_lp_solver_polyhedra(eltype(P); presolve=true))
    require(@__MODULE__, :Polyhedra; fun_name="chebyshev_center")
    if !ispolytopic(P)
        error("can only compute a Chebyshev center for polytopes")
    end

    Q = polyhedron(P; backend=backend)
    c, r = Polyhedra.chebyshevcenter(Q, solver)
    return c, r
end

function load_Polyhedra_polyhedron()
    return quote
        # see the interface file init_Polyhedra.jl for the imports

        """
            polyhedron(P::LazySet; [backend]=default_polyhedra_backend(P))

        Compute a set representation from `Polyhedra.jl`.

        ### Input

        - `P`       -- polyhedral set
        - `backend` -- (optional, default: call `default_polyhedra_backend(P)`)
                        the polyhedral computations backend

        ### Output

        A set representation in the `Polyhedra` library.

        ### Notes

        For further information on the supported backends see
        [Polyhedra's documentation](https://juliapolyhedra.github.io/).

        ### Algorithm

        This default implementation uses `tosimplehrep`, which computes the constraint
        representation of `P`. Set types preferring the vertex representation should
        implement their own method.
        """
        function polyhedron(P::LazySet; backend=default_polyhedra_backend(P))
            A, b = tosimplehrep(P)
            return Polyhedra.polyhedron(Polyhedra.hrep(A, b), backend)
        end
    end
end  # quote / load_Polyhedra_polyhedron()

function load_Polyhedra_GeometryBasics_triangulate_faces()
    return quote
        """
            triangulate_faces(X::LazySet)

        Triangulate the faces of a three-dimensional polytopic set.

        ### Input

        - `X` -- three-dimensional polytopic set

        ### Output

        A tuple `(p, c)` where `p` is a matrix, with each column containing a point, and
        `c` is a list of 3-tuples containing the indices of the points in each triangle.

        ### Notes

        This function triangulates all faces of a 3D polytope. The result is a list of (flat)
        triangles in 3D which describe the boundary of `X`.

        `X` must contain at least three vertices.
        """
        @validate function triangulate_faces(X::LazySet)
            P = polyhedron(X)
            mes = Mesh(P)
            coords = GeometryBasics.coordinates(mes)
            connection = GeometryBasics.faces(mes)

            ntriangles = length(connection)
            npoints = length(coords)
            @assert npoints == 3 * ntriangles
            points = Matrix{Float32}(undef, 3, npoints)

            for i in 1:npoints
                points[:, i] .= coords[i].data
            end

            connection_tup = getfield.(connection, :data)

            return points, connection_tup
        end
    end
end  # quote / load_Polyhedra_GeometryBasics_triangulate_faces()

"""
# Extended help

    isempty(P::LazySet, witness::Bool=false;
            [use_polyhedra_interface]::Bool=false, [solver]=nothing,
            [backend]=nothing)

### Input

- `witness` -- (optional, default: `false`) compute a witness if activated
- `use_polyhedra_interface` -- (optional, default: `false`) if `true`, we use
               the `Polyhedra` interface for the emptiness test
- `solver`  -- (optional, default: `nothing`) LP-solver backend; uses
               `default_lp_solver` if not provided
- `backend` -- (optional, default: `nothing`) backend for polyhedral
               computations in `Polyhedra`; uses `default_polyhedra_backend(P)`
               if not provided

### Notes

This implementation assumes that `P` is polyhedral.

The default value of the `backend` is set internally and depends on whether the
`use_polyhedra_interface` option is set.
If the option is set, we use `default_polyhedra_backend(P)`.

Witness production is not supported if `use_polyhedra_interface` is `true`.

### Algorithm

The algorithm sets up a feasibility LP for the constraints of `P`.
If `use_polyhedra_interface` is `true`, we call `Polyhedra.isempty`.
Otherwise, we set up the LP internally.
"""
function isempty(P::LazySet,
                 witness::Bool=false;
                 use_polyhedra_interface::Bool=false,
                 solver=nothing,
                 backend=nothing)
    @assert ispolyhedral(P) "this algorithm requires a polyhedral set"
    clist = constraints_list(P)
    if length(clist) < 2
        # catch corner case because of problems in LP solver for Rationals
        return witness ? (false, an_element(P)) : false
    end
    if use_polyhedra_interface
        return _isempty_polyhedron_polyhedra(P, witness; solver=solver,
                                             backend=backend)
    else
        return _isinfeasible(clist, witness; solver=solver)
    end
end

function _isempty_polyhedron(P::LazySet{N}, witness::Bool=false;
                             use_polyhedra_interface::Bool=false,
                             solver=nothing, backend=nothing) where {N}
    if use_polyhedra_interface
        return _isempty_polyhedron_polyhedra(P, witness; solver=solver,
                                             backend=backend)
    else
        return _isinfeasible(constraints_list(P), witness; solver=solver)
    end
end

function _isempty_polyhedron_polyhedra(P::LazySet{N}, witness::Bool=false;
                                       solver=nothing, backend=nothing) where {N}
    require(@__MODULE__, :Polyhedra; fun_name="isempty",
            explanation="with the active option `use_polyhedra_interface`")

    if isnothing(backend)
        backend = default_polyhedra_backend(P)
    end

    if isnothing(solver)
        result = Polyhedra.isempty(polyhedron(P; backend=backend))
    else
        result = Polyhedra.isempty(polyhedron(P; backend=backend), solver)
    end

    if result
        return _witness_result_empty(witness, true, N)
    elseif witness
        error("witness production is not supported yet")
    else
        return false
    end
end

"""
# Extended help

    linear_map(M::AbstractMatrix, P::LazySet; kwargs...)

### Algorithm

The default implementation assumes that `P` is polyhedral and applies an
algorithm based on the set type (see [`_linear_map_polyhedron`](@ref)).
"""
@validate function linear_map(M::AbstractMatrix, P::LazySet; kwargs...)
    if ispolyhedral(P)
        return _linear_map_polyhedron(M, P; kwargs...)
    else
        throw(ArgumentError("`linear_map` is not implemented for the given set"))
    end
end

function linear_map(αI::UniformScaling, X::LazySet)
    M = Diagonal(fill(αI.λ, dim(X)))
    return linear_map(M, X)
end

"""
# Extended help

    scale(α::Real, X::LazySet)

### Algorithm

The default implementation computes `linear_map` with the diagonal matrix ``α*I``.
"""
function scale(α::Real, X::LazySet)
    return linear_map(α * I, X)
end

function _scale_copy_inplace(α::Real, X::LazySet)
    Y = copy(X)
    scale!(α, Y)
    return Y
end

"""
    tohrep(P::LazySet)

Convert a polyhedral set to constraint representation.

### Input

- `P` -- polyhedral set

### Output

An `HPolyhedron` if `P` is bounded (via `isboundedtype`) or an `HPolytope`
otherwise.
"""
function tohrep(P::LazySet)
    @assert ispolyhedral(P) "cannot compute the constraint representation " *
                            "of non-polyhedral sets"

    T = isboundedtype(typeof(P)) ? HPolytope : HPolyhedron
    return T(constraints_list(P))
end

"""
    tovrep(P::LazySet)

Convert a polytopic set to vertex representation.

### Input

- `P` -- polytopic set

### Output

A `VPolytope`.
"""
function tovrep(P::LazySet)
    @assert ispolytopic(P) "cannot compute the vertex representation of non-polytopic sets"

    return VPolytope(vertices_list(P))
end

# The default implementation throws an error because Julia's default behavior
# leads to an error that is hard to understand.
@validate function ∈(::AbstractVector, X::LazySet)
    throw(ArgumentError("membership check for set type $(basetype(X)) is not implemented"))
end

function linear_map_inverse(A::AbstractMatrix, P::LazySet)
    @assert size(A, 1) == dim(P) "an inverse linear map of size $(size(A)) " *
                                 "cannot be applied to a set of dimension $(dim(P))"
    @assert ispolyhedral(P) "cannot compute the inverse linear map of " *
                            "non-polyhedral sets"
    return _affine_map_inverse(A, P)
end

function affine_map_inverse(A::AbstractMatrix, P::LazySet, b::AbstractVector)
    @assert size(A, 1) == dim(P) == length(b) "an inverse affine map of size $(size(A)) " *
                                              "and $(length(b)) cannot be applied to a " *
                                              "set of dimension $(dim(P))"
    @assert ispolyhedral(P) "cannot compute the inverse affine map of " *
                            "non-polyhedral sets"
    return _affine_map_inverse(A, P, b)
end

function _affine_map_inverse(A, P, b=nothing)
    constraints = _affine_map_inverse_hrep(A, P, b)
    if isempty(constraints)
        return Universe{eltype(P)}(size(A, 2))
    elseif length(constraints) == 2 && !isfeasible(constraints)
        return EmptySet{eltype(P)}(size(A, 2))
    end
    return HPolyhedron(constraints)
end

function vertices_list_1d(X::LazySet)
    l, h = extrema(X, 1)
    return _isapprox(l, h) ? [[l]] : [[l], [h]]
end

# internal function to detect parametric set types
isparametrictype(::Type{<:LazySet}) = false
