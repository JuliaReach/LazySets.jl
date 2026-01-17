# ========================
# Sampling from a LazySet
# ========================

"""
    AbstractSampler

Abstract type for defining new sampling methods.

### Notes

All subtypes should implement a `sample!(D, X, ::Method)` method where the
first argument is the output (vector of vectors), the second argument is the
set to be sampled, and the third argument is the sampler instance.
"""
abstract type AbstractSampler end

"""
    sample(X::LazySet{N}, num_samples::Int;
           [sampler]=_default_sampler(X),
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int, Nothing}=nothing,
           [include_vertices]=false,
           [VN]=Vector{N}) where {N}

Random sampling of an arbitrary set `X`.

### Input

- `X`           -- set to be sampled
- `num_samples` -- number of random samples
- `sampler`     -- (optional, default: `_default_sampler(X)`) the sampler used;
                   falls back to [`CombinedSampler`](@ref)
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding
- `include_vertices` -- (optional, default: `false`) option to include the
                   vertices of `X`
- `VN`          -- (optional, default: `Vector{N}`) vector type of the sampled
                   points

### Output

A vector of `num_samples` vectors.
If `num_samples` is not passed, the result is just one sample (not wrapped in a
vector).

### Algorithm

See the documentation of the respective `Sampler`.

### Notes

If `include_vertices == true`, we include all vertices computed with `vertices`.
Alternatively if a number ``k`` is passed, we plot the first ``k`` vertices
returned by `vertices(X)`.
"""
function sample(X::LazySet{N}, num_samples::Int;
                sampler=_default_sampler(X),
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int,Nothing}=nothing,
                include_vertices=false,
                VN=Vector{N}) where {N}
    D = Vector{VN}(undef, num_samples) # preallocate output
    sample!(D, X, sampler; rng=rng, seed=seed)

    if include_vertices != false
        k = (include_vertices isa Bool) ? Inf : include_vertices
        for v in vertices(X)
            push!(D, v)
            k -= 1
            if k <= 0
                break
            end
        end
    end

    return D
end

# without second argument: returns a single element instead of singleton list
function sample(X::LazySet; kwargs...)
    return sample(X, 1; kwargs...)[1]
end

# default sampling for LazySets
_default_sampler(::LazySet) = CombinedSampler()
function _default_sampler(::LineSegment{N}) where {N}
    return RejectionSampler(DefaultUniform(zero(N), one(N)), true, Inf, RejectionSamplerStatus())
end
_default_sampler(::HalfSpace) = HalfSpaceSampler()
_default_sampler(::Hyperplane) = HyperplaneSampler()
_default_sampler(::Line2D) = HyperplaneSampler()
_default_sampler(::AbstractSingleton) = SingletonSampler()
_default_sampler(::AbstractPolynomialZonotope) = PolynomialZonotopeSampler()
_default_sampler(::Universe) = UniverseSampler()
_default_sampler(H::AbstractHyperrectangle) = RejectionSampler(H)
_default_sampler(::AbstractZonotope) = PolynomialZonotopeSampler()

_rand(rng, U) = rand(rng, U)
_rand(rng, U::AbstractVector) = rand.(Ref(rng), U)

mutable struct RejectionSamplerStatus
    success::Bool

    function RejectionSamplerStatus()
        return new(false)
    end
end

"""
    RejectionSampler{D} <: AbstractSampler

Type used for rejection sampling of a bounded set `X`.

### Fields

- `distribution` -- (optional, default: `DefaultUniform`) distribution from
                    which the sample is drawn
- `tight`        -- (optional, default: `false`) set to `true` if the support of
                    the distribution is known to coincide with the set `X`
- `maxiter`      -- (optional, default: `Inf`) maximum number of iterations
                    before giving up

### Algorithm

Draw a sample ``x`` from a given distribution of a box-overapproximation of the
original set ``X`` in all ``n`` dimensions. The function rejects a drawn sample
``x`` and redraws as long as the sample is not contained in the original set
``X``, i.e., while ``x ∉ X``.

### Notes

The `maxiter` parameter is useful when sampling from sets that are small
compared to their box approximation, e.g., flat sets, for which the probability
of sampling from within the set is close to zero.
"""
struct RejectionSampler{D} <: AbstractSampler
    distribution::D
    tight::Bool
    maxiter::Number
    status::RejectionSamplerStatus
end

function get_success(S::RejectionSampler)
    return S.status.success
end

function set_success!(S::RejectionSampler)
    S.status.success = true
    return nothing
end

function RejectionSampler(distr; tight::Bool=false, maxiter=Inf)
    return RejectionSampler(distr, tight, maxiter, RejectionSamplerStatus())
end

function RejectionSampler(distr::DefaultUniform; tight::Bool=false, maxiter=Inf)
    return RejectionSampler([distr], tight, maxiter, RejectionSamplerStatus())
end

function RejectionSampler(X::LazySet, distribution=DefaultUniform;
                          tight::Bool=false, maxiter=Inf)
    # define the support of the distribution as the smallest box enclosing X
    n = dim(X)
    B = box_approximation(X)

    # distribution over B
    distr = [distribution(low(B, i), high(B, i)) for i in 1:n]

    return RejectionSampler(distr, tight, maxiter, RejectionSamplerStatus())
end

# the support of this distribution is always tight wrt. X
function RejectionSampler(X::AbstractHyperrectangle)
    n = dim(X)
    distr = [DefaultUniform(low(X, i), high(X, i)) for i in 1:n]
    return RejectionSampler(distr, true, Inf, RejectionSamplerStatus())
end

function sample!(D::Vector{VN}, X::LazySet, sampler::RejectionSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    U = sampler.distribution
    rng = reseed!(rng, seed)
    @inbounds for i in eachindex(D)
        w = _rand(rng, U)

        if !(sampler.tight)
            j = 1
            while w ∉ X && j <= sampler.maxiter
                w = _rand(rng, U)
                j += 1
            end
            if j > sampler.maxiter
                return D
            end
        end
        D[i] = w
    end
    set_success!(sampler)
    return D
end

function sample!(D::Vector{VN}, L::LineSegment,
                 sampler::RejectionSampler{<:DefaultUniform};
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    rng = reseed!(rng, seed)
    U = sampler.distribution
    @assert U.a >= zero(N) && U.b <= one(N) "illegal distribution domain"
    p = L.p
    q = L.q

    @inbounds for i in eachindex(D)
        λ = _rand(rng, U)
        D[i] = p + λ * (q - p)
    end
    set_success!(sampler)
    return D
end

"""
    RandomWalkSampler <: AbstractSampler

Type used for sampling from a convex polytope using its vertex representation.
This is especially useful if rejection sampling does not work because the
polytope is flat.

### Fields

- `variant` -- (optional, default: `true`) choice of a variant (see below)

### Notes

The sampling is not uniform - points in the center of the polytope are more
likely to be sampled.

The set to be sampled from must provide its vertices via `vertices_list`.

### Algorithm

Choose a random convex combination of the vertices of a convex polytope `X`.

If `variant == false`, we proceed as follows.
Let ``V = \\{v_i\\}_i`` denote the set of vertices of `X`.
Then any point ``p ∈ ℝ^n`` of the convex polytope ``X`` is a convex
combination of its vertices, i.e., ``p = ∑_{i} v_i α_i`` for some
(non-negative) coefficients ``\\{α_i\\}_i`` that add up to 1.
The algorithm chooses a random convex combination (the ``α_i``).
To produce this combination, we apply the finite-difference operator on a sorted
uniform sample over ``[0, 1]``; the method can be found in [Rubin81](@citet) and
[this StackExchange post](https://cs.stackexchange.com/a/3229).

If `variant == true`, we start from a random vertex and then repeatedly walk
toward a random vertex inside the polytope.
"""
struct RandomWalkSampler <: AbstractSampler
    variant::Bool
end

RandomWalkSampler() = RandomWalkSampler(true)

function sample!(D::Vector{VN}, X::LazySet, sampler::RandomWalkSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    rng = reseed!(rng, seed)
    U = DefaultUniform(zero(N), one(N))
    vlist = vertices_list(X)
    m = length(vlist)

    if sampler.variant
        @inbounds for i in eachindex(D)
            p = vlist[rand(1:m)]  # start from a random vertex
            for j in randperm(m)  # choose a random target vertex
                p += _rand(rng, U) * (vlist[j] - p)  # move toward next vertex
            end
            D[i] = p
        end
    else
        # vector used to store the combination coefficients
        r = Vector{N}(undef, m - 1)
        @inbounds for i in eachindex(D)
            # get a list of m uniform numbers (https://cs.stackexchange.com/a/3229)
            # and compute the corresponding linear combination in-place
            rand!(r, rng, U)
            sort!(r)
            D[i] = r[1] * vlist[1]  # r[1] - 0 == r[1]
            for j in 2:(m - 1)
                α = r[j] - r[j - 1]
                D[i] .+= α * vlist[j]
            end
            D[i] .+= (1 - r[m - 1]) * vlist[m]
        end
    end

    return D
end

"""
    CombinedSampler <: AbstractSampler

Type used for sampling arbitrary sets by trying different sampling strategies.

### Algorithm

The algorithm is to first try a [`RejectionSampler`](@ref) 10 times.
If that fails, it tries a [`RandomWalkSampler`](@ref).
"""
struct CombinedSampler <: AbstractSampler
    #
end

function sample!(D::Vector{VN}, X::LazySet, sampler::CombinedSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    # try rejection sampling 10 times
    tmp_sampler = RejectionSampler(X; maxiter=10)
    D2 = Vector{VN}(undef, 1)
    sample!(D2, X, tmp_sampler)
    if get_success(tmp_sampler)
        # it worked
        tmp_sampler = RejectionSampler(X)
        return sample!(D, X, tmp_sampler; rng=rng, seed=seed)
    end

    # try random-walk sampler
    tmp_sampler = RandomWalkSampler()
    try
        sample!(D, X, tmp_sampler; rng=rng, seed=seed)
        return D
    catch e
        @warn "sampling failed with the following error message:"
        rethrow(e)
    end
end

"""
    FaceSampler <: AbstractSampler

Type used for sampling from the `k`-faces of a set.

### Fields

- `dim` -- dimension of the faces to be sampled; a negative number is
           interpreted as `n - dim` where `n` is the dimension of the set

### Notes

For a three-dimensional polytope, the following face dimensions exist:
- 3-face – the polytope itself
- 2-faces – 2-dimensional polygonal faces
- 1-faces – 1-dimensional edges
- 0-faces – 0-dimensional vertices

For more information see
[Wikipedia](https://en.wikipedia.org/wiki/Face_(geometry)#k-face).

### Algorithm

Currently only hyperrectangles are supported.
For each point to be sampled, we randomly split the integers `1 .. n` into two
subgroups of size `k` and `n-k` respectively. For the i-th coordinate in the
first group, we sample in the interval `low(H, i) .. high(H, i)`. For the i-th
coordinate in the second group, we randomly pick either `low(H, i)` or
`high(H, i)`.
"""
struct FaceSampler <: AbstractSampler
    dim::Int
end

function sample!(D::Vector{VN}, X::LazySet, sampler::FaceSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    n = dim(X)
    k = sampler.dim

    # translate negative dimension
    if k < 0
        k += n
    end

    # full-dimensional sampling
    if k == n
        return sample!(D, X, _default_sampler(X); rng=rng, seed=seed)
    end

    0 <= k < n || throw(ArgumentError("cannot sample from " *
                                      "$(sampler.dim)-dimensional faces for a set of dimension $n"))

    rng = reseed!(rng, seed)
    return _sample_faces!(D, X, rng, k)
end

_choose_sorted(k, n; rng=GLOBAL_RNG) = sort!((1:n)[randperm(rng, n)][1:k])

function _sample_faces!(D::Vector{VN}, H::AbstractHyperrectangle, rng,
                        k) where {N,VN<:AbstractVector{N}}
    n = dim(H)

    @inbounds for j in eachindex(D)
        indices = _choose_sorted(k, n; rng=rng)
        # add one more index that is never visited
        indices = vcat(indices, 0)

        x = Vector{N}(undef, n)
        idx = 1
        for i in 1:n
            if indices[idx] == i
                # sample from interior in this dimension
                l = low(H, i)
                Δ = high(H, i) - l
                x[i] = Δ * rand(rng) + l
                idx += 1
            else
                # sample from border in this dimension
                x[i] = rand(rng, Bool) ? low(H, i) : high(H, i)
            end
        end
        D[j] = x
    end

    return D
end

"""
    HalfSpaceSampler{D} <: AbstractSampler

Type used for sampling from a half-space.

### Fields

- `distribution` -- (optional, default: `nothing`) distribution from which
                    samples are drawn

### Notes

If `distribution` is `nothing` (default), the sampling algorithm uses a
`DefaultUniform` over ``[0, 1]^n``.
"""
struct HalfSpaceSampler{D} <: AbstractSampler
    distribution::D
end

function HalfSpaceSampler()
    return HalfSpaceSampler(nothing)
end

function sample!(D::Vector{VN}, H::HalfSpace, sampler::HalfSpaceSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    rng = reseed!(rng, seed)
    U = sampler.distribution
    if isnothing(U)
        U = DefaultUniform(zero(N), one(N))
    end
    n = dim(H)
    @inbounds for i in eachindex(D)
        # sample a random point
        x = rand(rng, U, n)

        if x ∉ H
            # reflect point at the defining hyperplane
            x = _reflect_point_hyperplane(x, H.a, H.b)
        end
        D[i] = x
    end
end

"""
    HyperplaneSampler{D} <: AbstractSampler

Type used for sampling from a hyperplane.

### Fields

- `distribution` -- (optional, default: `nothing`) distribution from which
                    samples are drawn

### Notes

If `distribution` is `nothing` (default), the sampling algorithm uses a
`DefaultUniform` over ``[0, 1]^n``.
"""
struct HyperplaneSampler{D} <: AbstractSampler
    distribution::D
end

function HyperplaneSampler()
    return HyperplaneSampler(nothing)
end

function sample!(D::Vector{VN}, hp::Union{Hyperplane,Line2D},
                 sampler::HyperplaneSampler; rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    rng = reseed!(rng, seed)
    U = sampler.distribution
    if isnothing(U)
        U = DefaultUniform(zero(N), one(N))
    end
    n = dim(hp)
    @inbounds for i in eachindex(D)
        # sample a random point
        x = rand(rng, U, n)

        # project the point onto hp
        y = project(x, hp)
        D[i] = y
    end
end

"""
    SingletonSampler <: AbstractSampler

Type used for sampling from a singleton.
"""
struct SingletonSampler <: AbstractSampler
end

function sample!(D::Vector{VN}, S::AbstractSingleton, sampler::SingletonSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    x = element(S)
    @inbounds for i in eachindex(D)
        D[i] = copy(x)
    end
end

"""
    PolynomialZonotopeSampler{D} <: AbstractSampler

Type used for sampling from polynomial zonotopes.

### Fields

- `distribution` -- (optional, default: `nothing`) distribution from which
                    samples are drawn

### Notes

If `distribution` is `nothing` (default), the sampling algorithm uses a
`DefaultUniform` over ``[-1, 1]^n``.
"""
struct PolynomialZonotopeSampler{D} <: AbstractSampler
    distribution::D
end

function PolynomialZonotopeSampler()
    return PolynomialZonotopeSampler(nothing)
end

function sample!(D::Vector{VN}, P::Union{AbstractPolynomialZonotope,AbstractZonotope},
                 sampler::PolynomialZonotopeSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    rng = reseed!(rng, seed)
    U = sampler.distribution
    if isnothing(U)
        U = DefaultUniform(-one(N), one(N))
    end
    @inbounds for i in eachindex(D)
        # sample a random point
        x = copy(center(P))
        _add_generators!(x, P, U, rng)
        D[i] = x
    end
end

function _add_generators!(x, P::DensePolynomialZonotope, U, rng)
    p = size(P.E[1], 2)
    βs = rand(rng, U, p)  # dependent factors
    η = polynomial_order(P)
    for i in 1:η
        j = 1
        for coeffs in NondecreasingIndices(p, i)
            β = prod(βs[coeffs])
            if allequal(coeffs)
                # E
                Ei = P.E[i]
                jE = first(coeffs)  # column is any of the (equal) values in coeffs
                x .+= β * @view Ei[:, jE]
            else
                # F
                Fi = P.F[i - 1]  # offset of -1 because F's start with index 2
                x .+= β * @view Fi[:, j]
                j += 1
            end
        end
    end
    # G
    for g in eachcol(P.G)
        γ = _rand(rng, U)  # independent factors
        x .+= γ * g
    end
    return x
end

function _add_generators!(x, P::AbstractSparsePolynomialZonotope, U, rng)
    # dependent generators
    p = nparams(P)
    αₛ = rand(rng, U, p)
    γₛ = [prod(αₛ .^ e) for e in eachcol(expmat(P))]  # vector of monomial coefficients
    x .+= genmat_dep(P) * γₛ

    # independent generators
    q = ngens_indep(P)
    βₛ = rand(rng, U, q)
    x .+= genmat_indep(P) * βₛ

    return x
end

function _add_generators!(x, Z::Zonotope, U, rng)
    p = ngens(Z)
    ξ = rand(rng, U, p)
    x .+= genmat(Z) * ξ
    return x
end

"""
    UniverseSampler{D} <: AbstractSampler

Type used for sampling from universal sets.

### Fields

- `distribution` -- (optional, default: `nothing`) distribution from which
                    samples are drawn

### Notes

If `distribution` is `nothing` (default), the sampling algorithm uses a
`DefaultUniform` over ``[-100, 100]^n``.
"""
struct UniverseSampler{D} <: AbstractSampler
    distribution::D
end

function UniverseSampler()
    return UniverseSampler(nothing)
end

function sample!(D::Vector{VN}, X::LazySet, sampler::UniverseSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int,Nothing}=nothing) where {N,VN<:AbstractVector{N}}
    rng = reseed!(rng, seed)
    U = sampler.distribution
    if isnothing(U)
        U = DefaultUniform(N(-100), N(100))
    end
    n = dim(X)
    @inbounds for i in eachindex(D)
        # sample a random point
        D[i] = rand(rng, U, n)
    end
end

# =============================
# Code requiring Distributions
# =============================

"""
    _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                 [rng]::AbstractRNG=GLOBAL_RNG,
                                 [seed]::Union{Int, Nothing}=nothing) where {N}

Draw samples from a uniform distribution on an ``n``-dimensional unit sphere
using Muller's method.

### Input

- `D`    -- output, vector of points
- `n`    -- dimension of the sphere
- `p`    -- number of random samples
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

The modified vector `D`.

### Algorithm

This function implements Muller's method of normalized Gaussians [Muller59](@citet)
to uniformly sample over the ``n``-dimensional sphere ``S^n`` (which is the
bounding surface of the ``n``-dimensional unit ball).

Given ``n`` canonical Gaussian random variables ``Z₁, Z₂, …, Z_n``, the
distribution of the vectors

```math
\\dfrac{1}{α}\\left(z₁, z₂, …, z_n\\right)^T,
```
where ``α := \\sqrt{z₁² + z₂² + … + z_n²}``, is uniform over ``S^n``.
"""
function _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                      rng::AbstractRNG=GLOBAL_RNG,
                                      seed::Union{Int,Nothing}=nothing) where {N}
    return _sample_unit_nsphere_muller_distributions!(D, n, p; rng=rng, seed=seed)
end

function load_Distributions_sample()
    return quote
        using .Distributions: UnivariateDistribution

        function RejectionSampler(distr::UnivariateDistribution; tight::Bool=false)
            return RejectionSampler([distr]; tight=tight)
        end

        function _sample_unit_nsphere_muller_distributions!(D::Vector{Vector{N}},
                                                            n::Int, p::Int;
                                                            rng::AbstractRNG=GLOBAL_RNG,
                                                            seed::Union{Int,Nothing}=nothing) where {N}
            rng = reseed!(rng, seed)
            Zdims = [Distributions.Normal() for _ in 1:n] # normal distributions for each dimension
            v = Vector{N}(undef, n) # sample direction
            @inbounds for j in 1:p
                α = zero(N)
                for i in 1:n
                    v[i] = rand(rng, Zdims[i])
                    α += v[i]^2
                end
                D[j] = v ./ sqrt(α)
            end
            return D
        end
    end
end  # load_Distributions_sample
