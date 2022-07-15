import Base: rand

export Zonotope,
       scale,
       scale!,
       reduce_order,
       remove_zero_generators,
       remove_redundant_generators

using LazySets.Arrays: _vector_type, _matrix_type

"""
    Zonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}} <: AbstractZonotope{N}

Type that represents a zonotope.

### Fields

- `center`     -- center of the zonotope
- `generators` -- matrix; each column is a generator of the zonotope

### Notes

Mathematically, a zonotope is defined as the set

```math
Z = \\left\\{ x ∈ \\mathbb{R}^n : x = c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of *generators*.
This characterization defines a zonotope as the finite Minkowski sum of line
segments.
Zonotopes can be equivalently described as the image of a unit infinity-norm
ball in ``\\mathbb{R}^n`` by an affine transformation.

Zonotopes can be constructed in two different ways: either passing the generators
as a matrix, where each column represents a generator, or passing a list of vectors
where each vector represents a generator. Below we illustrate both ways.

### Examples

A two-dimensional zonotope with given center and set of generators:

```jldoctest zonotope_label
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> dim(Z)
2

julia> center(Z)
2-element Vector{Float64}:
 1.0
 0.0

julia> genmat(Z)
2×2 Matrix{Float64}:
 0.1  0.0
 0.0  0.1
```
Here, the first vector in the `Zonotope` constructor corresponds to the zonotope's
center, and each column of the second argument corresponds to a generator. The
functions `center` and `genmat` return the center and the generator matrix of this
zonotope respectively.

We can collect its vertices using `vertices_list`:

```jldoctest zonotope_label
julia> vertices_list(Z)
4-element Vector{Vector{Float64}}:
 [1.1, 0.1]
 [0.9, 0.1]
 [0.9, -0.1]
 [1.1, -0.1]
```

The support vector along a given direction can be computed using `σ`
(resp. the support function can be computed using `ρ`):

```jldoctest zonotope_label
julia> σ([1., 1.], Z)
2-element Vector{Float64}:
 1.1
 0.1
```

Zonotopes admit an alternative constructor that receives a list of
vectors, each vector representing a generator:

```jldoctest
julia> Z = Zonotope(ones(2), [[1., 0.], [0., 1.], [1., 1.]])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 1.0], [1.0 0.0 1.0; 0.0 1.0 1.0])

julia> genmat(Z)
2×3 Matrix{Float64}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
struct Zonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}} <: AbstractZonotope{N}
    center::VN
    generators::MN

    function Zonotope(center::VN, generators::MN) where {N,
                                                         VN<:AbstractVector{N},
                                                         MN<:AbstractMatrix{N}}
        @assert length(center) == size(generators, 1) "the dimension of the " *
            "center ($(length(center))) and the generators " *
            "($(size(generators, 1))) need to match"
        return new{N, VN, MN}(center, generators)
    end
end

isoperationtype(::Type{<:Zonotope}) = false
isconvextype(::Type{<:Zonotope}) = true

# constructor from center and list of generators
function Zonotope(center::VN, generators_list::AbstractVector{VN}) where {VN<:AbstractVector}
    G = to_matrix(generators_list, length(center))
    return Zonotope(center, G)
end

"""
    remove_zero_generators(Z::Zonotope)

Return a new zonotope removing the generators which are zero of the given zonotope.

### Input

- `Z` -- zonotope

### Output

If there are no zero generators, the result is the original zonotope `Z`.
Otherwise the result is a new zonotope that has the center and generators as `Z`
except for those generators that are zero.
"""
function remove_zero_generators(Z::Zonotope)
    G = Z.generators
    G2 = remove_zero_columns(G)
    if G === G2
        return Z
    end
    return Zonotope(Z.center, G2)
end

# --- AbstractCentrallySymmetric interface functions ---


"""
    center(Z::Zonotope)

Return the center of a zonotope.

### Input

- `Z` -- zonotope

### Output

The center of the zonotope.
"""
function center(Z::Zonotope)
    return Z.center
end


# --- AbstractZonotope interface functions ---


"""
    togrep(Z::Zonotope)

Return a generator representation of a zonotope.

### Input

- `Z` -- zonotope

### Output

The same set `Z`.
"""
function togrep(Z::Zonotope)
    return Z
end

"""
   genmat(Z::Zonotope)

Return the generator matrix of a zonotope.

### Input

- `Z` -- zonotope

### Output

A matrix where each column represents one generator of the zonotope `Z`.
"""
function genmat(Z::Zonotope)
    return Z.generators
end

"""
    generators(Z::Zonotope)

Return an iterator over the generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

An iterator over the generators of `Z`.
"""
function generators(Z::Zonotope)
    return generators_fallback(Z)
end

"""
    ngens(Z::Zonotope)

Return the number of generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

Integer representing the number of generators.
"""
ngens(Z::Zonotope) = size(Z.generators, 2)


# --- LazySet interface functions ---


"""
    rand(::Type{Zonotope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random zonotope.

### Input

- `Zonotope`       -- type for dispatch
- `N`              -- (optional, default: `Float64`) numeric type
- `dim`            -- (optional, default: 2) dimension
- `rng`            -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`           -- (optional, default: `nothing`) seed for reseeding
- `num_generators` -- (optional, default: `-1`) number of generators of the
                      zonotope (see comment below)

### Output

A random zonotope.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the argument `num_generators`.
For a negative value we choose a random number in the range `dim:2*dim` (except
if `dim == 1`, in which case we only create a single generator).
"""
function rand(::Type{Zonotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing,
              num_generators::Int=-1)
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    if num_generators < 0
        num_generators = (dim == 1) ? 1 : rand(dim:2*dim)
    end
    generators = randn(rng, N, dim, num_generators)
    return Zonotope(center, generators)
end


# --- Zonotope functions ---


"""
    scale(α::Real, Z::Zonotope)

Concrete scaling of a zonotope.

### Input

- `α` -- scalar
- `Z` -- zonotope

### Output

The zonotope obtained by applying the numerical scale to the center and
generators of ``Z``.
"""
function scale(α::Real, Z::Zonotope)
    c = α .* Z.center
    gi = α .* Z.generators
    return Zonotope(c, gi)
end

"""
    scale!(α::Real, Z::Zonotope)

Concrete scaling of a zonotope modifing `Z` in-place

### Input

- `α` -- scalar
- `Z` -- zonotope

### Output

The zonotope `Z` after applying the numerical scale `α` to its center and generators.
"""
function scale!(α::Real, Z::Zonotope)
    c = Z.center
    G = Z.generators
    c .= α .* c
    G .= α .* G
    return Z
end

# ==================================
# Zonotope order reduction methods
# ==================================

"""
    GIR05 <: AbstractReductionMethod

Zonotope order reduction method from [GIR05].

- [G05] A. Girard. *Reachability of Uncertain Linear Systems Using Zonotopes*, HSCC. Vol. 5. 2005.
"""
struct GIR05 <: AbstractReductionMethod end

"""
    COMB03 <: AbstractReductionMethod

Zonotope order reduction method from [COMB03].

- [COMB03] C. Combastel. *A state bounding observer based on zonotopes.* In Proc. of the European Control Conference, p. 2589–2594, 2003.
"""
struct COMB03 <: AbstractReductionMethod end

function reduce_order(Z::Zonotope, r::Number, method::AbstractReductionMethod=GIR05())
    r >= 1 || throw(ArgumentError("the target order should be at least 1, but it is $r"))
    c = Z.center
    G = Z.generators
    n, p = size(G)

    # r is bigger than the order of Z => don't reduce
    (r * n >= p) && return Z

    if isone(r)
        # if r = 1 => m = 0 and the generators need not be sorted
        Lred = _interval_hull(G, 1:p) 
        return Zonotope(c, Lred)
    end

    # sort generators
    indices = Vector{Int}(undef, p)
    _weighted_gens!(indices, G, method)

    # the first m generators have greatest weight
    m = floor(Int, n * (r - 1))

    # compute interval hull of L
    Lred = _interval_hull(G, view(indices, (m+1):p))

    # concatenate non-reduced and reduced generators
    Gred = _hcat_KLred(G, view(indices, 1:m), Lred)

    return Zonotope(c, Gred)
end

# Return the indices of the generators in G (= columns) sorted according to decreasing 2-norm.
# The generator index with highest score goes first.
function _weighted_gens!(indices, G::AbstractMatrix{N}, ::COMB03) where {N}
    p = size(G, 2)
    weights = Vector{N}(undef, p)
    @inbounds for j in 1:p
        v = view(G, :, j)
        weights[j] = norm(v, 2)
    end
    sortperm!(indices, weights, rev=true, initialized=false)
    return indices
end

# Return the indices of the generators in G (= columns) sorted according to ||⋅||₁ - ||⋅||∞ difference.
# The generator index with highest score goes first.
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
    sortperm!(indices, weights, rev=true, initialized=false)
    return indices
end

# compute interval hull of the generators of G (= columns) corresponding to `indices`
function _interval_hull(G::AbstractMatrix{N}, indices) where {N}
    n, p = size(G)
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
function _interval_hull(G::SMatrix{n, p, N, L}, indices) where {n, p, N, L}
    Lred = zeros(MMatrix{n, n, N})
    @inbounds for i in 1:n
        for j in indices
            Lred[i, i] += abs(G[i, j])
        end
    end
    return SMatrix{n, n}(Lred)
end

# implementation for static arrays
function _hcat_KLred(G::SMatrix{n, p, N, L1}, indices, Lred::SMatrix{n, n, N, L2}) where {n, p, N, L1, L2}
    m = length(indices)
    K = SMatrix{n, m}(view(G, :, indices))
    return hcat(K, Lred)
end

end end # quote / load_reduce_order_static

# ============================
# Zonotope splitting methods
# ============================

function _split(Z::Zonotope, j::Int)
    c, G = Z.center, Z.generators

    c₁ = similar(c)
    G₁ = similar(G)
    Z₁ = Zonotope(c₁, G₁)

    c₂ = similar(c)
    G₂ = similar(G)
    Z₂ = Zonotope(c₂, G₂)

    split!(Z₁, Z₂, Z, j)
end

function split!(Z₁::Zonotope, Z₂::Zonotope, Z::Zonotope, j::Int)
    c, G = Z.center, Z.generators
    n, p = size(G)
    @assert 1 <= j <= p "cannot split a zonotope with $p generators along index $j"

    c₁, G₁ = Z₁.center, Z₁.generators
    c₂, G₂ = Z₂.center, Z₂.generators
    copyto!(G₁, G)

    @inbounds for i in 1:n
        α = G[i, j] / 2
        c₁[i] = c[i] - α
        c₂[i] = c[i] + α
        G₁[i, j] = α
    end
    copyto!(G₂, G₁)

    return _split_ret(Z₁, Z₂)
end

_split_ret(Z₁::Zonotope, Z₂::Zonotope) = (Z₁, Z₂)

function load_split_static()
return quote

function _split_ret(Z₁::Zonotope{N, SV, SM}, Z₂::Zonotope{N, SV, SM}) where {N, n, p, SV<:MVector{n, N}, SM<:MMatrix{n, p, N}}
    Z₁ = Zonotope(SVector{n}(Z₁.center), SMatrix{n, p}(Z₁.generators))
    Z₂ = Zonotope(SVector{n}(Z₂.center), SMatrix{n, p}(Z₂.generators))
    return Z₁, Z₂
end

end end  # quote / load_split_static

function _split(Z::Zonotope, gens::AbstractVector, n::AbstractVector)
    p = length(gens)
    @assert p == length(n) "the number of generators doesn't match the " *
                            "number of indicated partitions ($p and $(length(n)))"

    @assert p <= ngens(Z) "the number of generators to split is greater " *
                          "than the number of generators of the zonotope ($p and $(ngens(Z)))"

    Zs = [Z]
    for (i, g) in enumerate(gens)
        for j = 1:n[i]
            km = length(Zs)
            for k = 1:km
                append!(Zs, split(Zs[k], g))
            end
            deleteat!(Zs, 1:km)
        end
    end
    return Zs
end

"""
    linear_map!(Zout::Zonotope, M::AbstractMatrix, Z::Zonotope)

Compute the concrete linear map of a zonotope storing the result in `Zout`.

### Input

- `Zout` -- zonotope (output)
- `M`    -- matrix
- `Z`    -- zonotope

### Output

The zonotope `Zout`, which is modified in-place.
"""
function linear_map!(Zout::Zonotope, M::AbstractMatrix, Z::Zonotope)
    mul!(Zout.center, M, Z.center)
    mul!(Zout.generators, M, Z.generators)
    return Zout
end

"""
    _bound_intersect_2D(Z::Zonotope, L::Line2D)

Return the support function in the direction [0, 1] of the intersection between
the given zonotope and line.

### Input

- `Z` -- zonotope
- `L` -- vertical line 2D

### Output

The support function in the direction [0, 1] of the intersection between the
given zonotope and line.

### Notes

The algorithm assumes that the given line is vertical and that the intersection
between the given sets is not empty.

### Algorithm

This function implements [Algorithm 8.2, 1].

[1] *Colas Le Guernic. Reachability Analysis of Hybrid Systems with Linear
Continuous Dynamics. Computer Science [cs]. Université Joseph-Fourier - Grenoble
I, 2009. English. fftel-00422569v2f*
"""
function _bound_intersect_2D(Z::Zonotope, L::Line2D)
    c = center(Z)
    P = copy(c)
    G = genmat(Z)
    r = ngens(Z)
    g(x) = view(G, :, x)
    for i = 1:r
        gi = g(i)
        if !_isupwards(gi)
            gi .= -gi
        end
        P .= P - gi
    end
    G = sortslices(G, dims=2, by=x->atan(x[2], x[1])) # sort gens
    if P[1] < L.b
        G .= G[:,end:-1:1]
    end
    j = 1
    while isdisjoint(LineSegment(P, P+2g(j)), L)
        P .= P + 2g(j)
        j += 1
        if j > size(G, 2)
            error("Got unexpected error, check that the sets intersect")
        end
    end
    singleton = intersection(LineSegment(P, P+2g(j)), L)
    return element(singleton)[2]
end
# ====================================
# Zonotope vertex enumeration methods
# ====================================

function _vertices_list_2D(c::AbstractVector{N}, G::AbstractMatrix{N}; apply_convex_hull::Bool) where {N}
    if same_sign(G)
        return _vertices_list_2D_positive(c, G, apply_convex_hull=apply_convex_hull)
    else
        # FIXME generalized 2D vertices list function is not implemented yet
        # See LazySets#2209
        return _vertices_list_iterative(c, G, apply_convex_hull=apply_convex_hull)
    end
end

function _vertices_list_2D_positive(c::AbstractVector{N}, G::AbstractMatrix{N}; apply_convex_hull::Bool) where {N}
    n, p = size(G)

    # TODO special case p = 1 or p = 2 ?

    sorted_G = sortslices(G, dims=2, by=x->atan(x[2], x[1]))
    index = ones(N, p, 2*p)
    @inbounds for i in 1:p
        index[i, i+1:i+p-1] .= -one(N)
    end
    index[:, 1] .= -one(N)
    V = sorted_G * index .+ c
    vlist = [V[:, i] for i in 1:2*p]

    if apply_convex_hull
        convex_hull!(vlist)
    end
    return vlist
end

function _vertices_list_iterative(c::VN, G::MN; apply_convex_hull::Bool) where {N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}}
    p = size(G, 2)
    vlist = Vector{VN}()
    sizehint!(vlist, 2^p)

    for ξi in Iterators.product([(1, -1) for i = 1:p]...)
        push!(vlist, c .+ G * collect(ξi))
    end

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

# special case 2D zonotope of order 1/2
function _vertices_list_2D_order_one_half(c::VN, G::MN; apply_convex_hull::Bool) where {N, VN<:AbstractVector{N}, MN}
    vlist = Vector{VN}(undef, 2)
    g = view(G, :, 1)
    @inbounds begin
        vlist[1] = c .+ g
        vlist[2] = c .- g
    end
    return apply_convex_hull ? _two_points_2d!(vlist) : vlist
end

# special case 2D zonotope of order 1
function _vertices_list_2D_order_one(c::VN, G::MN; apply_convex_hull::Bool) where {N, VN<:AbstractVector{N}, MN}
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

"""
    remove_redundant_generators(Z::Zonotope{N}) where {N}

Remove all redundant (pairwise linearly dependent) generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

A new zonotope with fewer generators, or the same zonotope if no generator could
be removed.

### Algorithm

For each generator ``g_j`` that has not been checked yet, we find all other
generators that are linearly dependent with ``g_j``.
Then we combine those generators into a single generator.

For one-dimensional zonotopes we use a more efficient implementation where we
just take the absolute sum of all generators.
"""
function remove_redundant_generators(Z::Zonotope{N}) where {N}
    if dim(Z) == 1  # more efficient implementation in 1D
        return _remove_redundant_generators_1d(Z)
    end

    G = genmat(Z)
    G = remove_zero_columns(G)
    p = size(G, 2)
    removed_zero_generators = p < ngens(Z)
    deleted = false
    done = falses(p)
    G_new = _vector_type(typeof(G))[]  # list of new column vectors
    @inbounds for j1 in 1:p
        if done[j1]  # skip if the generator was already removed
            continue
        end
        # "done[j1] = true" not needed because we will never look at it again
        gj1 = G[:, j1]
        for j2 in (j1+1):p  # look at all generators to the right
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
        G_new = reduce(hcat, G_new)  # convert list of column vectors to matrix
        return Zonotope(center(Z), G_new)
    elseif removed_zero_generators
        return Zonotope(center(Z), G)
    end
    return Z  # return the original zonotope if no generator was removed
end

function _remove_redundant_generators_1d(Z)
    G = genmat(Z)
    g = sum(abs, G)
    return Zonotope(center(Z), hcat(g))
end

"""
    low(Z::Zonotope, i::Int)

Return the lower coordinate of a zonotope in a given dimension.

### Input

- `Z` -- zonotope
- `i` -- dimension of interest

### Output

The lower coordinate of the zonotope in the given dimension.
"""
function low(Z::Zonotope, i::Int)
    G = genmat(Z)
    v = center(Z, i)
    @inbounds for j in 1:ngens(Z)
        v -= abs(G[i, j])
    end
    return v
end

"""
    high(Z::Zonotope, i::Int)

Return the higher coordinate of a zonotope in a given dimension.

### Input

- `Z` -- zonotope
- `i` -- dimension of interest

### Output

The higher coordinate of the zonotope in the given dimension.
"""
function high(Z::Zonotope, i::Int)
    G = genmat(Z)
    v = center(Z, i)
    @inbounds for j in 1:ngens(Z)
        v += abs(G[i, j])
    end
    return v
end
