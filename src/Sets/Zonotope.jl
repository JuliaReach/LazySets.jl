import Base: rand,
             split

export Zonotope,
       scale,
       reduce_order

"""
    Zonotope{N<:Real} <: AbstractZonotope{N}

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

Zonotopes can be constructed in two different ways: either passing the generators as a matrix, where
each column represents a generator, or passing a list of vectors where each vector represents a generator.
Below we illustrate both ways.

The optional argument `remove_zero_generators` controls whether we remove zero
columns from the `generators` matrix. This option is active by default.

### Examples

A two-dimensional zonotope with given center and set of generators:

```jldoctest zonotope_label
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> dim(Z)
2
```
Here, each column of the second input corresponds to a generator.

We can collect its vertices using `vertices_list`:

```jldoctest zonotope_label
julia> vertices_list(Z)
4-element Array{Array{Float64,1},1}:
 [1.1, 0.1]
 [0.9, 0.1]
 [0.9, -0.1]
 [1.1, -0.1]
```

The support vector along a given direction can be computed using `σ`
(resp. the support function can be computed using `ρ`):

```jldoctest zonotope_label
julia> σ([1., 1.], Z)
2-element Array{Float64,1}:
 1.1
 0.1
```

Zonotopes admit an alternative constructor that receives a list of
vectors, each vector representing a generator:

```jldoctest
julia> Z = Zonotope(ones(2), [[1., 0.], [0., 1.], [1., 1.]])
Zonotope{Float64}([1.0, 1.0], [1.0 0.0 1.0; 0.0 1.0 1.0])

julia> Z.generators
2×3 Array{Float64,2}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
struct Zonotope{N<:Real} <: AbstractZonotope{N}
    center::AbstractVector{N}
    generators::AbstractMatrix{N}

    function Zonotope(center::AbstractVector{N}, generators::AbstractMatrix{N};
                      remove_zero_generators::Bool=true) where {N<:Real}
        @assert length(center) == size(generators, 1) "the dimension of the " *
            "center ($(length(center))) and the generators " *
            "($(size(generators, 1))) need to match"
        if remove_zero_generators
            generators = delete_zero_columns!(generators)
        end
        new{N}(center, generators)
    end
end

isoperationtype(::Type{<:Zonotope}) = false

# constructor from center and list of generators
Zonotope(center::AbstractVector{N}, generators_list::AbstractVector{VN};
         remove_zero_generators::Bool=true
        ) where {N<:Real, VN<:AbstractVector{N}} =
    Zonotope(center, hcat(generators_list...);
             remove_zero_generators=remove_zero_generators)


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(Z::Zonotope{N})::Vector{N} where {N<:Real}

Return the center of a zonotope.

### Input

- `Z` -- zonotope

### Output

The center of the zonotope.
"""
function center(Z::Zonotope{N})::Vector{N} where {N<:Real}
    return Z.center
end


# --- AbstractZonotope interface functions ---


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
    ngens(Z::Zonotope)::Int

Return the number of generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

Integer representing the number of generators.
"""
ngens(Z::Zonotope)::Int = size(Z.generators, 2)


# --- LazySet interface functions ---


"""
    rand(::Type{Zonotope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Zonotope{N}

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
              num_generators::Int=-1
             )::Zonotope{N}
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
    reduce_order(Z::Zonotope, r)::Zonotope

Reduce the order of a zonotope by overapproximating with a zonotope with less
generators.

### Input

- `Z` -- zonotope
- `r` -- desired order

### Output

A new zonotope with less generators, if possible.

### Algorithm

This function implements the algorithm described in A. Girard's
*Reachability of Uncertain Linear Systems Using Zonotopes*, HSCC. Vol. 5. 2005.

If the desired order is smaller than one, the zonotope is *not* reduced.
"""
function reduce_order(Z::Zonotope{N}, r)::Zonotope{N} where {N<:Real}
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
    split(Z::Zonotope, j::Int)

Return two zonotopes obtained by splitting the given zonotope.

### Input

- `Z` -- zonotope
- `j` -- index of the generator to be split

### Output

The zonotope obtained by splitting `Z` into two zonotopes such that
their union is `Z` and their intersection is possibly non-empty.

### Algorithm

This function implements [Prop. 3, 1], that we state next. The zonotope
``Z = ⟨c, g^{(1, …, p)}⟩`` is split into:

```math
Z₁ = ⟨c - \\frac{1}{2}g^{(j)}, (g^{(1, …,j-1)}, \\frac{1}{2}g^{(j)}, g^{(j+1, …, p)})⟩ \\\\
Z₂ = ⟨c + \\frac{1}{2}g^{(j)}, (g^{(1, …,j-1)}, \\frac{1}{2}g^{(j)}, g^{(j+1, …, p)})⟩,
```
such that ``Z₁ ∪ Z₂ = Z`` and ``Z₁ ∩ Z₂ = Z^*``, where

```math
Z^* = ⟨c, (g^{(1,…,j-1)}, g^{(j+1,…, p)})⟩.
```

[1] *Althoff, M., Stursberg, O., & Buss, M. (2008). Reachability analysis of
nonlinear systems with uncertain parameters using conservative linearization.
In Proc. of the 47th IEEE Conference on Decision and Control.*
"""
function split(Z::Zonotope, j::Int)
    @assert 1 <= j <= ngens(Z) "cannot split a zonotope with $(ngens(Z)) generators along index $j"
    c, G = Z.center, Z.generators
    Gj = G[:, j]
    Gj_half = Gj / 2

    c₁ = c - Gj_half
    c₂ = c + Gj_half

    G₁ = copy(G)
    G₁[:, j] = Gj_half
    G₂ = copy(G₁)

    Z₁ = Zonotope(c₁, G₁)
    Z₂ = Zonotope(c₂, G₂)
    return Z₁, Z₂
end
