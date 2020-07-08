import Base: rand,
             split

export Zonotope,
       scale,
       reduce_order,
       remove_zero_generators

using LazySets.Arrays: _vector_type, _matrix_type
using StaticArrays: SMatrix, SVector, MMatrix, MVector

"""
    Zonotope{N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}} <: AbstractZonotope{N}

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
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> dim(Z)
2

julia> center(Z)
2-element Array{Float64,1}:
 1.0
 0.0

julia> genmat(Z)
2×2 Array{Float64,2}:
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
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.0, 1.0], [1.0 0.0 1.0; 0.0 1.0 1.0])

julia> genmat(Z)
2×3 Array{Float64,2}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
struct Zonotope{N<:Real, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}} <: AbstractZonotope{N}
    center::VN
    generators::MN

    function Zonotope(center::VN, generators::MN) where {N<:Real,
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
function Zonotope(center::VN, generators_list::AbstractVector{VN}) where {N<:Real, VN<:AbstractVector{N}}
    MT = _matrix_type(VN)
    G = MT(undef, length(center), length(generators_list))
    for (j, gj) in enumerate(generators_list)
        @inbounds G[:, j] = gj
    end
    return Zonotope(center, G)
end

"""
    remove_zero_generators(Z::Zonotope{N, VN, MN}) where {N<:Real,
                                                          VN<:AbstractVector{N},
                                                          MN<:AbstractMatrix{N}}

Return a new zonotope removing the generators which are zero of the given zonotope.

### Input

- `Z` -- zonotope

### Output

If there are no zero generators, the result is the original zonotope `Z`.
Otherwise the result is a new zonotope that has the center and generators as `Z`
except for those generators that are zero.
"""
function remove_zero_generators(Z::Zonotope{N, VN, MN}) where {N<:Real,
                                                               VN<:AbstractVector{N},
                                                               MN<:AbstractMatrix{N}}
    G = Z.generators
    G2 = remove_zero_columns(G)
    if G === G2
        return Z
    end
    return Zonotope(Z.center, G2)
end

# --- AbstractCentrallySymmetric interface functions ---


"""
    center(Z::Zonotope{N}) where {N<:Real}

Return the center of a zonotope.

### Input

- `Z` -- zonotope

### Output

The center of the zonotope.
"""
function center(Z::Zonotope{N}) where {N<:Real}
    return Z.center
end


# --- AbstractZonotope interface functions ---


"""
    togrep(Z::Zonotope)

Return a generator representation of a zonotope.

### Input

- `Z` -- zonotopic set

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
    reduce_order(Z::Zonotope, r::Union{Integer, Rational})

Reduce the order of a zonotope by overapproximating with a zonotope with less
generators.

### Input

- `Z` -- zonotope
- `r` -- desired order

### Output

A new zonotope with less generators, if possible.

### Algorithm

See `overapproximate(Z::Zonotope{N}, ::Type{<:Zonotope}, r::Union{Integer, Rational}) where {N<:Real}` for details.
"""
function reduce_order(Z::Zonotope{N}, r::Union{Integer, Rational}) where {N<:Real}
    return overapproximate(Z, Zonotope, r)
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
    copyto!(G₂, G)

    @inbounds for i in 1:n
        α = G[i, j] / 2
        c₁[i] = c[i] - α
        c₂[i] = c[i] + α
        G₁[i, j] = α
        G₂[i, j] = α
    end
    return Z₁, Z₂
end

"""
    split(Z::Zonotope{N, SVector{n, N}, <:SMatrix{n, p, N}}, j::Int) where {N, n, p}

Return two zonotopes obtained by splitting the given zonotope.

### Input

- `Z` -- zonotope
- `j` -- index of the generator to be split

### Output

The zonotopes obtained by splitting `Z` into two zonotopes such that
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
function split(Z::Zonotope{N, SVector{n, N}, <:SMatrix{n, p, N}}, j::Int) where {N, n, p}
    @assert 1 <= j <= p "cannot split a zonotope with $p generators along index $j"
    c, G = Z.center, Z.generators

    c₁ = MVector{n, N}(undef)
    c₂ = MVector{n, N}(undef)

    G₁ = MMatrix{n, p}(G)
    G₂ = MMatrix{n, p}(G)

    @inbounds for i in 1:n
        α = G[i, j] / 2
        c₁[i] = c[i] - α
        c₂[i] = c[i] + α
        G₁[i, j] = α
        G₂[i, j] = α
    end

    Z₁ = Zonotope(SVector{n}(c₁), SMatrix{n, p}(G₁))
    Z₂ = Zonotope(SVector{n}(c₂), SMatrix{n, p}(G₂))
    return Z₁, Z₂
end

"""
    split(Z::Zonotope, gens::AbstractVector{Int}, n::AbstractVector{Int})

Return a vector zonotopes obtained form splitting in the given generators the
given zonotope.

### Input

- `Z` -- zonotope
- `gens` -- vector of indices of the generators to be splitted
- `n` -- vector of integers describing the number of partitions in the
         corresponding generator

### Output

The zonotopes obtained by splitting `Z` into `2^{n_i}` zonotopes for each
generator `i` such that their union is `Z` and their intersection is
possibly non-empty.
"""
function split(Z::Zonotope, gens::AbstractVector, n::AbstractVector)
    @assert length(gens) == length(n) "the number of generators doesn't match the" *
    " number of indicated partitions ($(length(gens)) and $(length(n)))"
    @assert length(gens) <= ngens(Z) "the number of generators to split is greater" *
    " than the number of generators of the zonotope (($(length(gens)) and $(ngens(Z)))"
    Zs = [Z]
    for i = 1:length(gens)
        for j = 1:n[i]
            km = length(Zs)
            for k = 1:km
                append!(Zs, split(Zs[k], gens[i]))
            end
            deleteat!(Zs, 1:km)
        end
    end
    return Zs
end
