import Base: split

export AbstractZonotope,
       genmat,
       generators,
       togrep,
       reduce_order,
       remove_redundant_generators

"""
    AbstractZonotope{N} <: AbstractCentrallySymmetricPolytope{N}

Abstract type for zonotopic sets.

### Notes

Mathematically, a zonotope is defined as the set

```math
Z = \\left\\{ c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i ∈ [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c ∈ ℝ^n`` is its center and ``\\{g_i\\}_{i=1}^p``,
``g_i ∈ ℝ^n``, is the set of generators.
This characterization defines a zonotope as the finite Minkowski sum of line
segments.
Zonotopes can be equivalently described as the image of a unit infinity-norm
ball in ``ℝ^n`` by an affine transformation.

See [`Zonotope`](@ref) for a standard implementation of this interface.

Every concrete `AbstractZonotope` must define the following functions:

- `generators(::AbstractZonotope)` -- return an iterator over the generators
- `genmat(::AbstractZonotope)` -- return a generator matrix

Since the functions `genmat` and `generators` can be defined in terms of each
other, it is sufficient to only genuinely implement one of them and let the
implementation of the other function call the fallback implementation
`genmat_fallback` resp. `generators_fallback`.

The subtypes of `AbstractZonotope` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractZonotope)
5-element Vector{Any}:
 AbstractHyperrectangle
 HParallelotope
 LineSegment
 Zonotope
 ZonotopeMD
```
"""
abstract type AbstractZonotope{N} <: AbstractCentrallySymmetricPolytope{N} end

"""
    generators(Z::AbstractZonotope)

Return an iterator over the generators of a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

An iterator over the generators of `Z`.
"""
function generators(::AbstractZonotope) end

"""
    genmat(Z::AbstractZonotope)

Return a generator matrix of a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

A generator matrix of `Z`.
"""
function genmat(::AbstractZonotope) end

"""
    genmat_fallback(Z::AbstractZonotope; [gens]=generators(Z), [ngens]=nothing)

Fallback definition of `genmat` for zonotopic sets.

### Input

- `Z`     -- zonotopic set
- `gens`  -- (optional; default: `generators(Z)`) iterator over generators
- `ngens` -- (optional; default: `nothing`) number of generators or `nothing` if
             unknown

### Output

A matrix where each column represents one generator of `Z`.

### Notes

Passing the number of generators (`ngens`) is more efficient, since otherwise
the generators have to be obtained from the iterator (`gens`) and stored in an
intermediate vector until the final result matrix can be allocated.
"""
function genmat_fallback(Z::AbstractZonotope; gens=generators(Z), ngens=nothing)
    if isempty(gens)
        N = eltype(Z)
        return Matrix{N}(undef, dim(Z), 0)
    elseif isnothing(ngens)
        return _genmat_fallback_generic(Z, gens)
    else
        return _genmat_fallback_ngens(Z, gens, ngens)
    end
end

function _genmat_fallback_generic(Z::AbstractZonotope, gens)
    gens = collect(gens)
    return _genmat_fallback_ngens(Z, gens, length(gens))
end

function _genmat_fallback_ngens(Z::AbstractZonotope{N}, gens, ngens) where {N}
    G = Matrix{N}(undef, dim(Z), ngens)
    for (i, g) in enumerate(gens)
        G[:, i] = g
    end
    return G
end

"""
    generators_fallback(Z::AbstractZonotope)

Fallback definition of `generators` for zonotopic sets.

### Input

- `Z` -- zonotopic set

### Output

An iterator over the generators of `Z`.
"""
function generators_fallback(Z::AbstractZonotope)
    return eachcol(genmat(Z))
end

"""
    togrep(Z::AbstractZonotope)

Return a generator representation of a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

The same set in generator representation.
This fallback implementation returns a `Zonotope`; however, more specific
implementations may return other generator representations.
"""
function togrep(Z::AbstractZonotope)
    return convert(Zonotope, Z)
end

"""
    ngens(Z::AbstractZonotope)

Return the number of generators of a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

An integer representing the number of generators.
"""
function ngens(Z::AbstractZonotope)
    return length(generators(Z))
end

"""
    order(Z::AbstractZonotope)

Return the order of a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

A rational number representing the order of the zonotopic set.

### Notes

The order of a zonotopic set is defined as the quotient of its number of
generators and its dimension.
"""
function order(Z::AbstractZonotope)
    return ngens(Z) // dim(Z)
end

"""
# Extended help

    ρ(d::AbstractVector, Z::AbstractZonotope)

### Algorithm

The support value is ``cᵀ d + ‖Gᵀ d‖₁``, where ``c`` is the center and ``G`` is
the generator matrix of `Z`.
"""
@validate function ρ(d::AbstractVector, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    return dot(c, d) + abs_sum(d, G)
end

"""
# Extended help

    σ(d::AbstractVector, Z::AbstractZonotope)

### Notes

If the direction has norm zero, the vertex with ``ξ_i = 1 \\ \\ ∀ i = 1,…, p``
is returned.
"""
@validate function σ(d::AbstractVector, Z::AbstractZonotope)
    G = genmat(Z)
    v = G * sign_cadlag.(At_mul_B(G, d))
    v .+= center(Z)
    return v
end

"""
# Extended help

    in(x::AbstractVector, Z::AbstractZonotope; solver=nothing)

### Input

- `solver` -- (optional, default: `nothing`) the backend used to solve the
              linear program

### Examples

```jldoctest
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1]);

julia> [1.0, 0.2] ∈ Z
false
julia> [1.0, 0.1] ∈ Z
true
```

### Notes

If `solver == nothing`, we fall back to `default_lp_solver(N)`.

### Algorithm

The membership problem is reduced to the following linear program.
Let ``p`` and ``n`` be the number of generators and ambient dimension,
respectively.
We consider the ``p``-dimensional space of elements ``(ξ_1, …, ξ_p)``
constrained to ``ξ_i ∈ [-1, 1]`` for all ``i = 1, …, p`` such that
``x-c = Gξ`` holds.
"""
@validate function in(x::AbstractVector, Z::AbstractZonotope; solver=nothing)
    res, _ = _in_AbstractZonotope_LP(x, Z; solver=solver)
    return res
end

# result: `(res::Bool, lp)` where `lp` can be either the LP solver output or `nothing`
function _in_AbstractZonotope_LP(x::AbstractVector, Z::AbstractZonotope; solver=nothing)
    p = ngens(Z)
    if p == 0
        # no generators can cause trouble in LP solver
        return (_isapprox(x, center(Z)), nothing)
    end

    N = promote_type(eltype(x), eltype(Z))
    A = genmat(Z)
    sense = '='
    b = x .- center(Z)
    lbounds = -one(N)
    ubounds = one(N)
    obj = zeros(N, p)

    if isnothing(solver)
        solver = default_lp_solver(N)
    end
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    res = is_lp_optimal(lp.status)  # Infeasible or Unbounded => false
    return (res, lp)
end

"""
# Extended help

    linear_map(M::AbstractMatrix, Z::AbstractZonotope)

### Output

A `Zonotope`.

### Algorithm

We apply the linear map to the center and the generators.

If the map has output dimension 1, a specialized algorithm ensures that the
resulting zonotope only has a single generator (or none if the map is zero).
"""
@validate function linear_map(M::AbstractMatrix, Z::AbstractZonotope)
    if size(M, 1) == 1
        # yields only one generator
        return _linear_map_zonotope_1D(M, Z)
    else
        return _linear_map_zonotope_nD(M, Z)
    end
end

function _linear_map_zonotope_1D(M::AbstractMatrix, Z::LazySet)
    N = promote_type(eltype(M), eltype(Z))
    c = M * center(Z)
    gi = zero(N)
    @inbounds for g in generators(Z)
        gi += abs(sum(M[1, i] * g[i] for i in eachindex(g)))
    end
    G = iszero(gi) ? Matrix{N}(undef, 1, 0) : hcat(gi)
    return Zonotope(c, G)
end

function _linear_map_zonotope_nD(M::AbstractMatrix, Z::LazySet)
    c = M * center(Z)
    G = M * genmat(Z)
    return Zonotope(c, G)
end

"""
# Extended help

    vertices_list(Z::AbstractZonotope; [apply_convex_hull]::Bool=true)

### Input

- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         computation with the convex hull of the points

### Algorithm

#### Two-dimensional case

We use a trick to speed up enumerating vertices of 2-dimensional zonotopic
sets with all generators in the first quadrant or third quadrant (same sign).
Namely, sort the generators by angle and add them clockwise in increasing
order and counterclockwise in decreasing order. A more detailed explanation can
be found [here](https://math.stackexchange.com/q/3356460).

To avoid the cumulative sum from both directions separately, we build a 2D index
matrix to sum generators for both directions in one matrix-vector product.

#### General case

If the zonotopic set has ``p`` generators, each vertex is the result of summing
the center with some linear combination of generators, where the combination
factors are ``ξ_i ∈ \\{-1, 1\\}``.

There are at most ``2^p`` distinct vertices. Use the flag `apply_convex_hull` to
control whether a convex-hull algorithm is applied to the vertices computed by
this method; otherwise, redundant vertices may be present.
"""
function vertices_list(Z::AbstractZonotope; apply_convex_hull::Bool=true)
    c = center(Z)
    G = genmat(Z)
    n, p = size(G)

    # empty generators => sole vertex is the center
    if p == 0
        return [c]
    end

    if n == 1
        return vertices_list(convert(Interval, Z))

    elseif n == 2
        if p == 1
            return _vertices_list_zonotope_2D_order_one_half(c, G;
                                                             apply_convex_hull=apply_convex_hull)
        elseif p == 2
            return _vertices_list_zonotope_2D_order_one(c, G; apply_convex_hull=apply_convex_hull)
        else
            return _vertices_list_zonotope_2D(c, G; apply_convex_hull=apply_convex_hull)
        end

    else
        Gred = remove_zero_columns(G)
        return _vertices_list_zonotope_iterative(c, Gred; apply_convex_hull=apply_convex_hull)
    end
end

function _vertices_list_zonotope_2D(c::AbstractVector{N}, G::AbstractMatrix{N};
                                    apply_convex_hull::Bool) where {N}
    if same_sign(G)
        if any(e -> e < 0, G)
            G = -G
        end
        return _vertices_list_zonotope_2D_positive(c, G; apply_convex_hull=apply_convex_hull)
    else
        # TODO generalized 2D vertices list function is not implemented yet
        # See LazySets#2209
        return _vertices_list_zonotope_iterative(c, G; apply_convex_hull=apply_convex_hull)
    end
end

function _vertices_list_zonotope_2D_positive(c::AbstractVector{N}, G::AbstractMatrix{N};
                                             apply_convex_hull::Bool) where {N}
    n, p = size(G)

    sorted_G = sortslices(G; dims=2, by=x -> atan(x[2], x[1]))
    index = ones(N, p, 2 * p)
    @inbounds for i in 1:p
        index[i, (i + 1):(i + p - 1)] .= -one(N)
    end
    index[:, 1] .= -one(N)
    V = sorted_G * index .+ c
    vlist = [V[:, i] for i in 1:(2 * p)]

    if apply_convex_hull
        convex_hull!(vlist)
    end
    return vlist
end

function _vertices_list_zonotope_iterative(c::VN, G::AbstractMatrix{N};
                                           apply_convex_hull::Bool) where {N,VN<:AbstractVector{N}}
    p = size(G, 2)
    vlist = Vector{VN}()
    sizehint!(vlist, 2^p)

    for ξi in Iterators.product([(1, -1) for i in 1:p]...)
        push!(vlist, c .+ G * collect(ξi))
    end

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

# special case 2D zonotope of order 1/2
function _vertices_list_zonotope_2D_order_one_half(c::VN, G::AbstractMatrix{N};
                                                   apply_convex_hull::Bool) where {N,
                                                                                   VN<:AbstractVector{N}}
    vlist = Vector{VN}(undef, 2)
    g = view(G, :, 1)
    @inbounds begin
        vlist[1] = c .+ g
        vlist[2] = c .- g
    end
    return apply_convex_hull ? _two_points_2d!(vlist) : vlist
end

# special case 2D zonotope of order 1
function _vertices_list_zonotope_2D_order_one(c::VN, G::AbstractMatrix{N};
                                              apply_convex_hull::Bool) where {N,
                                                                              VN<:AbstractVector{N}}
    vlist = Vector{VN}(undef, 4)
    a = [one(N), one(N)]
    b = [one(N), -one(N)]
    @inbounds begin
        vlist[1] = c .+ G * a
        vlist[2] = c .- G * a
        vlist[3] = c .+ G * b
        vlist[4] = c .- G * b
    end
    return apply_convex_hull ? _four_points_2d!(vlist) : vlist
end

# TODO this function is currently not used
function _vertices_list_2D(c::AbstractVector{N}, G::AbstractMatrix{N};
                           apply_convex_hull::Bool) where {N}
    if apply_convex_hull
        return _vertices_list_zonotope_iterative(c, G; apply_convex_hull=apply_convex_hull)
    end

    angles = mapslices(_angles, G; dims=1)[1, :]
    perm = sortperm(angles)
    sorted_angles = angles[perm]
    sorted_G = G[:, perm]
    polygons = Vector{VPolygon{N}}()
    sizehint!(polygons, 4)

    @inbounds for i in zip(0:90:360, 90:90:360)
        index = i[1] .<= sorted_angles .< i[2]
        if sum(index) > 0
            push!(polygons, _single_quadrant_vertices_enum(sorted_G[:, index], true))
        end
    end

    return vertices_list(translate(reduce(minkowski_sum, polygons), c))
end

@inline function _angles(point::AbstractVector)
    return (atand(point[2], point[1]) + 360) % 360
end

function _single_quadrant_vertices_enum(G::AbstractMatrix, sorted::Bool=false)
    if !sorted
        G = sortslices(G; dims=2, by=_angles)
    end
    return VPolygon(2 * cumsum(hcat(G, -G); dims=2) .- sum(G; dims=2))
end

"""
# Extended help

    constraints_list(P::AbstractZonotope)

### Algorithm

This is the (inefficient) fallback implementation for rational numbers.
It first computes the vertices and then converts the corresponding polytope
to constraint representation.
"""
function constraints_list(Z::AbstractZonotope)
    return _constraints_list_vrep(Z)
end

@inline function _constraints_list_vrep(Z::AbstractZonotope)
    return constraints_list(VPolytope(vertices_list(Z)))
end

"""
# Extended help

    constraints_list(Z::AbstractZonotope{<:AbstractFloat})

### Notes

The main algorithm assumes that the generator matrix is full rank.
The result has ``2 \\binom{p}{n-1}`` (with ``p`` being the number of generators
and ``n`` being the ambient dimension) constraints, which is optimal under this
assumption.
If this assumption is not given, the implementation tries to work around.

### Algorithm

We follow the algorithm presented in [AlthoffSB10](@citet). Three cases are not covered by that
algorithm, so we handle them separately. The first case is zonotopes in one
dimension. The second case is that there are fewer generators than dimensions,
``p < n``, or the generator matrix is not full rank, in which case we fall back
to the (slower) computation based on the vertex representation. The third case
is that the zonotope is flat in some dimensions, in which case we project the
zonotope to the non-flat dimensions and extend the result later.
"""
function constraints_list(Z::AbstractZonotope{<:AbstractFloat})
    return _constraints_list_zonotope(Z)
end

function _constraints_list_zonotope(Z::AbstractZonotope{N}) where {N<:AbstractFloat}
    n = dim(Z)

    # special handling of the 1D case
    if n == 1
        return _constraints_list_1d(Z)
    end

    p = ngens(Z)

    # check whether to use the fallback implementation in V-rep
    use_vrep = (p < n)  # order < 1
    if !use_vrep
        # remove redundant generators and check again
        Z = remove_redundant_generators(Z)
        p = ngens(Z)
        use_vrep = (p < n)
    end
    if use_vrep
        return _constraints_list_vrep(Z)
    end

    extend_indices = flat_dimensions(Z)
    isflat_Z = !isempty(extend_indices)
    if isflat_Z
        Z_orig = Z
        Z = project(Z, setdiff(1:n, extend_indices))
        n = dim(Z)
    end

    G = genmat(Z)
    Gᵀ = transpose(G)
    c = center(Z)
    m = binomial(p, n - 1)
    constraints = Vector{HalfSpace{N,Vector{N}}}()
    sizehint!(constraints, 2m)
    for columns in StrictlyIncreasingIndices(p, n - 1)
        c⁺ = cross_product(view(G, :, columns))
        iszero(c⁺) && continue
        normalize!(c⁺, 2)

        Δd = sum(abs, Gᵀ * c⁺)

        c⁺c = dot(c⁺, c)
        d⁺ = c⁺c + Δd
        c⁻ = -c⁺
        d⁻ = -c⁺c + Δd

        if isflat_Z
            c⁺ = extend_with_zeros(c⁺, extend_indices)
            c⁻ = extend_with_zeros(c⁻, extend_indices)
        end

        push!(constraints, HalfSpace(c⁺, d⁺))
        push!(constraints, HalfSpace(c⁻, d⁻))
    end
    if isflat_Z
        # add constraints about flat dimensions
        c = center(Z_orig)
        n = length(c)
        @inbounds for i in extend_indices
            ci = c[i]
            a = zeros(N, n)
            a[i] = one(N)
            push!(constraints, HalfSpace(a, ci))
            a = zeros(N, n)
            a[i] = -one(N)
            push!(constraints, HalfSpace(a, -ci))
        end
    end
    return constraints
end

function _constraints_list_1d(Z::AbstractZonotope{N}) where {N}
    c = center(Z, 1)
    g = sum(abs, genmat(Z))
    return [HalfSpace([N(1)], c + g), HalfSpace([N(-1)], g - c)]
end

function flat_dimensions(Z::AbstractZonotope)
    G = genmat(Z)
    @static if VERSION >= v"1.7"
        res = findall(iszero, eachrow(G))
    else
        res = findall(iszero, collect(eachrow(G)))
    end
    return res
end

"""
    split(Z::AbstractZonotope, j::Int)

Return two zonotopes obtained by splitting the given zonotopic set.

### Input

- `Z` -- zonotopic set
- `j` -- index of the generator to be split

### Output

The zonotope obtained by splitting `Z` into two zonotopes such that
their union is `Z` and their intersection is possibly non-empty.

### Algorithm

This function implements [AlthoffSB08; Prop. 3](@citet), which we state next.
The zonotopic set ``Z = ⟨c, g^{(1, …, p)}⟩`` is split into:

```math
Z₁ = ⟨c - \\frac{1}{2}g^{(j)}, (g^{(1, …,j-1)}, \\frac{1}{2}g^{(j)}, g^{(j+1, …, p)})⟩ \\\\
Z₂ = ⟨c + \\frac{1}{2}g^{(j)}, (g^{(1, …,j-1)}, \\frac{1}{2}g^{(j)}, g^{(j+1, …, p)})⟩,
```
such that ``Z₁ ∪ Z₂ = Z`` and ``Z₁ ∩ Z₂ = Z^*``, where

```math
Z^* = ⟨c, (g^{(1,…,j-1)}, g^{(j+1,…, p)})⟩.
```
"""
function split(Z::AbstractZonotope, j::Int)
    return _split(convert(Zonotope, Z), j)
end

"""
    split(Z::AbstractZonotope, gens::AbstractVector{Int},
          nparts::AbstractVector{Int})

Split a zonotopic set along the given generators into a vector of zonotopes.

### Input

- `Z`    -- zonotopic set
- `gens` -- vector of indices of the generators to be split
- `n`    -- vector of integers describing the number of partitions in the
            corresponding generator

### Output

The zonotopes obtained by splitting `Z` into `2^{n_i}` zonotopes for each
generator `i` such that their union is `Z` and their intersection is
possibly non-empty.

### Examples

Splitting of a two-dimensional zonotopic set along its first generator:

```jldoctest zonotope_label
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> split(Z, [1], [1])
2-element Vector{Zonotope{Float64, Vector{Float64}, Matrix{Float64}}}:
 Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.95, 0.0], [0.05 0.0; 0.0 0.1])
 Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.05, 0.0], [0.05 0.0; 0.0 0.1])
```
Here, the first vector in the arguments corresponds to the zonotopic set's
generator to be split, and the second vector corresponds to the exponent of
`2^n` parts that the set will be split into along the corresponding generator.

As an example, below we split a two-dimensional zonotope along both of its
generators, each time into four parts.

```
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> split(Z, [1, 2], [2, 2])
16-element Vector{Zonotope{Float64, Vector{Float64}, Matrix{Float64}}}:
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.925, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.925, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.925, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.925, 0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.975, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.975, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.975, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.975, 0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.025, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.025, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.025, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.025, 0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.075, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.075, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.075, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.075, 0.075], [0.025 0.0; 0.0 0.025])
```
"""
function split(Z::AbstractZonotope, gens::AbstractVector{Int},
               nparts::AbstractVector{Int})
    return _split(convert(Zonotope, Z), gens, nparts)
end

# project via the concrete linear map, which is typically efficient
@validate function project(Z::AbstractZonotope{N}, block::AbstractVector{Int}; kwargs...) where {N}
    n = dim(Z)
    M = projection_matrix(block, n, N)
    Z2 = linear_map(M, Z)

    if get(kwargs, :remove_zero_generators, true)
        Z2 = remove_zero_generators(Z2)
    end
    return Z2
end

"""
    remove_redundant_generators(Z::AbstractZonotope)

Remove all redundant (pairwise linearly dependent) generators of a zonotopic
set.

### Input

- `Z` -- zonotopic set

### Output

A new zonotope with fewer generators, or the same zonotopic set if no generator
could be removed.

### Algorithm

By default this implementation returns the input zonotopic set. Subtypes of
`AbstractZonotope` whose generators can be removed have to define a new method.
"""
function remove_redundant_generators(Z::AbstractZonotope)
    return Z  # fallback implementation
end

function remove_redundant_generators(G::AbstractMatrix{N}) where {N<:Rational}
    if size(G, 1) == 1  # more efficient implementation in 1D
        return _remove_redundant_generators_1d(G)
    end

    G = remove_zero_columns(G)
    p = size(G, 2)
    deleted = false
    done = falses(p)
    G_new = vector_type(typeof(G))[]  # list of new column vectors
    @inbounds for j1 in 1:p
        if done[j1]  # skip if the generator was already removed
            continue
        end
        # "done[j1] = true" not needed because we will never look at it again
        gj1 = G[:, j1]
        for j2 in (j1 + 1):p  # look at all generators to the right
            if done[j2]  # skip if the generator was already removed
                continue
            end
            gj2 = G[:, j2]
            answer, factor = ismultiple(gj1, gj2)
            if answer
                # column j2 is a multiple of column j1
                if factor > zero(N)
                    gj1 += gj2
                else
                    gj1 -= gj2
                end
                done[j2] = true
                deleted = true
            end
        end
        push!(G_new, gj1)
    end

    if deleted
        return reduce(hcat, G_new)  # convert list of column vectors to matrix
    end
    return G
end

function remove_redundant_generators(G::AbstractMatrix)
    if size(G, 1) == 1  # more efficient implementation in 1D
        return _remove_redundant_generators_1d(G)
    end

    G = remove_zero_columns(G)
    p = size(G, 2)

    if p <= 1
        return G
    end

    # normalize columns
    norms = Vector{eltype(G)}(undef, p)
    Gnorm = similar(G)
    @inbounds for (j, col) in enumerate(eachcol(G))
        norms[j] = norm(col)
        unit_col = col ./ norms[j]
        idx = findfirst(!iszero, unit_col)

        # flip the vector if the first element is negative
        if !isnothing(idx) && unit_col[idx] < 0
            Gnorm[:, j] = -unit_col
        else
            Gnorm[:, j] = unit_col
        end
    end

    # sort column in ascending order
    cols = eachcol(Gnorm)
    if VERSION < v"1.9"
        cols = collect(cols)
    end
    ord = sortperm(cols)
    merged = Vector{Vector{eltype(G)}}()

    this = ord[1]
    cur_dir = view(Gnorm, :, this)
    cur_sum = norms[this]  # accumulated norm of identical cols

    # traverse the sorted matrix and compare cols
    @inbounds for j in 2:p
        prev = this
        this = ord[j]

        if _isapprox(view(Gnorm, :, prev), view(Gnorm, :, this))
            cur_sum += norms[this]
        else
            push!(merged, cur_dir .* cur_sum) # merge multiples

            cur_dir = view(Gnorm, :, this)
            cur_sum = norms[this]
        end
    end

    push!(merged, cur_dir .* cur_sum)

    G_new = hcat(merged...)
    return G_new
end

function _remove_redundant_generators_1d(G::AbstractMatrix)
    g = sum(abs, G)
    return hcat(g)
end

"""
    AbstractReductionMethod

Abstract supertype for order-reduction methods of a zonotopic set.
"""
abstract type AbstractReductionMethod end

"""
    GIR05 <: AbstractReductionMethod

Zonotope order-reduction method from [Girard05](@citet).
"""
struct GIR05 <: AbstractReductionMethod end

"""
    COMB03 <: AbstractReductionMethod

Zonotope order-reduction method from [Combastel03](@citet).
"""
struct COMB03 <: AbstractReductionMethod end

"""
    ASB10 <: AbstractReductionMethod

Zonotope order-reduction method from [AlthoffSB10](@citet).
"""
struct ASB10 <: AbstractReductionMethod end

"""
    SRMB16 <: AbstractReductionMethod

Zonotope order-reduction method from [ScottRMB16](@citet).

### Fields

- `ϵ` -- (optional; default: `1e-6`) pivot threshold
- `δ` -- (optional; default: `1e-3`) volume threshold

### Notes

The method reorders the generator matrix using reduced row echelon form (rref) to the form
``[T ~ V]``, then iteratively removes one generator from ``V`` while updating ``T``.
"""
struct SRMB16{N<:Number} <: AbstractReductionMethod
    ϵ::N
    δ::N

    SRMB16(ϵ::N=1e-6, δ::N=1e-3) where {N<:Number} = new{N}(ϵ, δ)
end

"""
    reduce_order(Z::AbstractZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())

Reduce the order of a zonotopic set by overapproximating with a zonotope with
fewer generators.

### Input

- `Z`      -- zonotopic set
- `r`      -- desired order
- `method` -- (optional, default: `GIR05()`) the reduction method used

### Output

A new zonotope with fewer generators, if possible.

### Algorithm

The available algorithms are:

```jldoctest; setup = :(using LazySets: subtypes, AbstractReductionMethod)
julia> subtypes(AbstractReductionMethod)
4-element Vector{Any}:
 LazySets.ASB10
 LazySets.COMB03
 LazySets.GIR05
 LazySets.SRMB16
```

See the documentation of each algorithm for references. Most methods split the
given zonotopic set `Z` into two zonotopes, `K` and `L`, where `K` contains the
most "representative" generators and `L` contains the generators that are
reduced, `Lred`, using a box overapproximation. This methodology varies slightly
for [`SRMB16`](@ref). We follow the notation from [YangS18](@citet). See also
[KopetzkiSA17](@citet).
"""
function reduce_order(Z::AbstractZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    r >= 1 || throw(ArgumentError("the target order should be at least 1, " *
                                  "but it is $r"))
    n = dim(Z)
    p = ngens(Z)

    # if r is bigger than the order of Z => do not reduce
    (r * n >= p) && return Z
    return _reduce_order_zonotope_common(Z, r, n, p, method)
end

function _reduce_order_zonotope_common(Z, r, n, p, method::Union{ASB10,COMB03,GIR05})
    c = center(Z)
    G = genmat(Z)

    if isone(r)
        # if r = 1 => m = 0 and the generators need not be sorted
        Lred = _approximate_reduce_order(c, G, 1:p, method)
        return Zonotope(c, Lred)
    end

    # sort generators
    indices = Vector{Int}(undef, p)
    _weighted_gens!(indices, G, method)

    # the first m generators have greatest weight
    m = floor(Int, n * (r - 1))

    return _absorb_generators_zonotope(c, G, indices, m, method)
end

# absorb all but the first `m` generators (according to the sorting in `indices`) with a box
function _absorb_generators_zonotope(c, G, indices, m, method)

    # compute interval hull of L
    Lred = _approximate_reduce_order(c, G, view(indices, (m + 1):length(indices)), method)

    # concatenate non-reduced and reduced generators
    Gred = _hcat_KLred(G, view(indices, 1:m), Lred)

    return Zonotope(c, Gred)
end

function _reduce_order_zonotope_common(Z, r, n, p, method::SRMB16)
    c = center(Z)
    G = genmat(Z)

    # reorder columns: G ↦ [T V] and also obtain IR = [I R] where R = T⁻¹V
    # TODO if G is rank deficient, fall back to COMB03
    TV, IR = _factorG(G, method.ϵ, method.δ)
    T = TV[:, 1:n]
    V = TV[:, (n + 1):end]
    R = IR[:, (n + 1):end]

    # iteratively remove one generator v by overapproximating with a parallelotope;
    m = floor(Int, n * r)
    while p > m
        # choose v that minimizes the volume error of this approximation
        j_min = argmin([prod(1 .+ abs.(R[:, j])) - sum(abs, R[:, j]) for j in axes(R, 2)])

        # slice out column j_min
        V = V[:, [1:(j_min - 1); (j_min + 1):end]]
        R₋ = R[:, [1:(j_min - 1); (j_min + 1):end]]

        diag_matrix = Diagonal(1 .+ abs.(R[:, j_min]))
        T = T * diag_matrix
        if !isempty(R)
            R = inv(diag_matrix) * R₋
        end

        p -= 1
    end

    G_red = hcat(T, V)
    return Zonotope(c, G_red)
end

# reorder columns of G as [T V] such that T is invertible; also return IR of the
# same shape and such that IR = [I R] where R = T⁻¹V
function _factorG(G::AbstractMatrix, ϵ::N, δ::N) where {N}
    TV = copy(G)
    IR = copy(G)
    n, ng = size(G)

    # rref with column swap
    for k in 1:n
        # normalize rows k:n
        for i in k:n
            row_norm = norm(IR[i, :], 1)
            if row_norm > ϵ
                IR[i, :] ./= row_norm
            end
        end

        submatrix = IR[k:n, k:ng]
        max_val, idx = findmax(abs.(submatrix))
        i_max = k + idx[1] - 1
        j_max = k + idx[2] - 1

        if abs(max_val) ≤ ϵ
            break
        end

        # swap rows and columns
        if i_max != k
            IR[[k, i_max], :] = IR[[i_max, k], :]
        end

        if j_max != k
            TV[:, [k, j_max]] = TV[:, [j_max, k]]
            IR[:, [k, j_max]] = IR[:, [j_max, k]]
        end

        _rref!(IR)
    end

    #  extra column swaps until all |IR_ij| ≤ 1 + δ
    while true
        (max_val, max_idx) = findmax(abs.(IR))
        (max_val ≤ 1 + δ) && break
        k, j = Tuple(max_idx)

        TV[:, [k, j]] = TV[:, [j, k]]
        IR[:, [k, j]] = IR[:, [j, k]]

        pivot = IR[k, k]
        for i in 1:n
            i == k && continue
            factor = IR[i, k] / pivot
            IR[i, :] -= factor .* IR[k, :]
        end
        IR[k, :] ./= pivot
    end

    return TV, IR
end

# reduced row echelon form
function _rref!(A::AbstractMatrix)
    nr, nc = size(A)
    pivot = 1

    for r in 1:nr
        if pivot > nc
            return A
        end

        i = r
        while i ≤ nr && A[i, pivot] == 0
            i += 1
            if i > nr
                i = r
                pivot += 1
                if pivot > nc
                    return A
                end
            end
        end

        A[[i, r], :] = A[[r, i], :]

        if A[r, pivot] != 0
            A[r, :] ./= A[r, pivot]
        end

        for i in 1:nr
            i == r && continue
            A[i, :] -= A[i, pivot] * A[r, :]
        end

        pivot += 1
    end

    return A
end

# approximate with a box
function _approximate_reduce_order(c, G, indices, ::Union{COMB03,GIR05})
    return _interval_hull(G, indices)
end

# approximate with a parallelotope
function _approximate_reduce_order(c, G, indices, ::ASB10)
    Ztilde = Zonotope(c, view(G, :, indices))
    Ψ = Approximations._overapproximate_hparallelotope(Ztilde)
    return genmat(Ψ)
end

# Return the indices of the generators in G (= columns) sorted according to
# decreasing 2-norm. The generator index with highest score goes first.
function _weighted_gens!(indices, G::AbstractMatrix{N}, ::Union{ASB10,COMB03}) where {N}
    p = size(G, 2)
    weights = Vector{N}(undef, p)
    @inbounds for j in 1:p
        v = view(G, :, j)
        weights[j] = norm(v, 2)
    end
    sortperm!(indices, weights; rev=true, initialized=false)
    return indices
end

# Return the indices of the generators in G (= columns) sorted according to
# ‖⋅‖₁ - ‖⋅‖∞ difference. The generator index with highest score goes first.
function _weighted_gens!(indices, G::AbstractMatrix{N}, ::GIR05) where {N}
    n, p = size(G)
    weights = Vector{N}(undef, p)
    @inbounds for j in 1:p
        aux_norm_1 = zero(N)
        aux_norm_inf = zero(N)
        for i in 1:n
            abs_Gij = abs(G[i, j])
            aux_norm_1 += abs_Gij
            if aux_norm_inf < abs_Gij
                aux_norm_inf = abs_Gij
            end
        end
        weights[j] = aux_norm_1 - aux_norm_inf
    end
    sortperm!(indices, weights; rev=true, initialized=false)
    return indices
end

# compute interval hull of the generators of G (= columns) corresponding to
# `indices`
function _interval_hull(G::AbstractMatrix{N}, indices) where {N}
    n = size(G, 1)
    Lred = zeros(N, n, n)
    @inbounds for i in 1:n
        for j in indices
            Lred[i, i] += abs(G[i, j])
        end
    end
    return Lred
end

# given an n x p matrix G and a vector of m integer indices with m <= p,
# concatenate the columns of G given by `indices` with the matrix Lred
function _hcat_KLred(G::AbstractMatrix, indices, Lred::AbstractMatrix)
    K = view(G, :, indices)
    return hcat(K, Lred)
end

function load_reduce_order_static()
    return quote

        # implementation for static arrays
        function _interval_hull(G::SMatrix{n,p,N,L}, indices) where {n,p,N,L}
            Lred = zeros(MMatrix{n,n,N})
            @inbounds for i in 1:n
                for j in indices
                    Lred[i, i] += abs(G[i, j])
                end
            end
            return SMatrix{n,n}(Lred)
        end

        # implementation for static arrays
        function _hcat_KLred(G::SMatrix{n,p,N,L1}, indices,
                             Lred::SMatrix{n,n,N,L2}) where {n,p,N,L1,L2}
            m = length(indices)
            K = SMatrix{n,m}(view(G, :, indices))
            return hcat(K, Lred)
        end
    end
end # quote / load_reduce_order_static

"""
# Extended help

    reflect(Z::AbstractZonotope)

### Algorithm

If ``Z`` has center ``c`` and generator matrix ``G``, then ``-Z`` has center
``-c`` and generator matrix ``G``. For the latter, observe that ``G`` and ``-G``
behave the same way.
"""
function reflect(Z::AbstractZonotope)
    return Zonotope(-center(Z), genmat(Z))
end

# internal function; defined here due to dependency StaticArrays and submodules
function _genmat_static(::AbstractZonotope) end

"""
    _norm_1(Z::AbstractZonotope)

Compute the exact ``ℓ₁`` norm of a zonotope with generator matrix ``G ∈ \\mathbb{R}^{d×n}``

### Notes

The implementation exploits that the mapping ``ξ ↦ ‖ c + ∑_{i=1}^n ξ_i g_i ‖_1`` is
a convex function of the coefficients ``ξ_i``.  As a result, its maximum over the hypercube
``[-1,1]^n`` is attained at one of the ``2^n`` corners.

This algorithm has time complexity ``\\mathcal{O}(2ⁿ · d)``, and thus is only practical for
zonotopes with a small number of generators.
"""
function _norm_1(Z::AbstractZonotope)
    N = eltype(Z)
    c = center(Z)
    G = genmat(Z)
    n = size(G, 2)
    dirs = DiagDirections(n)
    norm = N(-Inf)

    @inbounds for v in dirs
        aux = sum(abs.((c + G * v)))
        norm = max(norm, aux)
    end

    return norm
end

function _norm_Inf(Z::AbstractZonotope)
    c = center(Z)
    r = _box_radius(Z)
    res = zero(eltype(Z))
    @inbounds for (ci, ri) in zip(c, r)
        # dominant box edge in dimension i
        if ci > 0
            vi = abs(ci + ri)
        else
            vi = abs(ci - ri)
        end
        res = max(res, vi)
    end
    return res
end

function linear_map_inverse(A::AbstractMatrix, Z::AbstractZonotope)
    @assert size(A, 1) == dim(Z) "an inverse linear map of size $(size(A)) " *
                                 "cannot be applied to a set of dimension $(dim(Z))"
    return _affine_map_inverse_zonotope(A, Z)
end

function affine_map_inverse(A::AbstractMatrix, Z::AbstractZonotope, b::AbstractVector)
    @assert size(A, 1) == dim(Z) == length(b) "an inverse affine map of size $(size(A)) " *
                                              "and $(length(b)) cannot be applied to a " *
                                              "set of dimension $(dim(Z))"
    return _affine_map_inverse_zonotope(A, Z, b)
end

# Z = {y = c + GB}
# A(c₀ + G₀B) + b = c + GB
# if A is invertible:
# <=>  c₀ = A⁻¹(c - b) ∧ G₀ = A⁻¹G
function _affine_map_inverse_zonotope(A::AbstractMatrix, Z::AbstractZonotope,
                                      b::Union{AbstractVector,Nothing}=nothing)
    if !isinvertible(A)
        # result may not be representable as a zonotope -> use default implementation for polyhedra
        return _affine_map_inverse(A, Z, b)
    end
    Ainv = inv(A)
    c = Ainv * center(Z)
    if !isnothing(b)
        c .-= Ainv * b
    end
    G = Ainv * genmat(Z)
    return Zonotope(c, G)
end

function _box_radius(Z::AbstractZonotope)
    return sum(abs, genmat(Z); dims=2)[:]
end

function high(Z::AbstractZonotope)
    return center(Z) .+ _box_radius(Z)
end

function high(Z::AbstractZonotope, i::Int)
    return center(Z, i) .+ _box_radius(Z)[i]
end

function low(Z::AbstractZonotope)
    return center(Z) .- _box_radius(Z)
end

function low(Z::AbstractZonotope, i::Int)
    return center(Z, i) .- _box_radius(Z)[i]
end
