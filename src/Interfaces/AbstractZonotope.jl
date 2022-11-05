import Base: ∈, split

export AbstractZonotope,
       genmat,
       generators,
       ngens,
       order,
       translate,
       translate!,
       togrep,
       split!,
       reduce_order

"""
    AbstractZonotope{N} <: AbstractCentrallySymmetricPolytope{N}

Abstract type for zonotopic sets.

### Notes

Mathematically, a zonotope is defined as the set

```math
Z = \\left\\{ c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its center and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of generators.
This characterization defines a zonotope as the finite Minkowski sum of line
segments.
Zonotopes can be equivalently described as the image of a unit infinity-norm
ball in ``\\mathbb{R}^n`` by an affine transformation.

See [`Zonotope`](@ref) for a standard implementation of this interface.

Every concrete `AbstractZonotope` must define the following functions:
- `genmat(::AbstractZonotope{N})` -- return the generator matrix
- `generators(::AbstractZonotope{N})` -- return an iterator over the generators

Since the functions `genmat` and `generators` can be defined in terms of each
other, it is sufficient to only genuinely implement one of them and let the
implementation of the other function call the fallback implementation
`genmat_fallback` resp. `generators_fallback`.

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractZonotope)
5-element Vector{Any}:
 AbstractHyperrectangle
 HParallelotope
 LineSegment
 RotatedHyperrectangle
 Zonotope
```
"""
abstract type AbstractZonotope{N} <: AbstractCentrallySymmetricPolytope{N} end

isconvextype(::Type{<:AbstractZonotope}) = true

# --- fallback implementations of AbstractZonotope functions ---


"""
    genmat_fallback(Z::AbstractZonotope{N};
                    [gens]=generators(Z),
                    [ngens]=nothing) where {N}

Fallback definition of `genmat` for zonotopic sets.

### Input

- `Z`     -- zonotopic set
- `gens`  -- (optional; default: `generators(Z)`) iterator over generators
- `ngens` -- (optional; default: `nothing`) number of generators or `nothing` if
             unknown

### Output

A matrix where each column represents one generator of `Z`.

### Notes

Passing the number of generators is much more efficient as otherwise the
generators have to be obtained from the iterator (`gens`) and stored in an
intermediate vector until the final result matrix can be allocated.
"""
function genmat_fallback(Z::AbstractZonotope{N};
                         gens=generators(Z),
                         ngens=nothing) where {N}
    if isempty(gens)
        return Matrix{N}(undef, dim(Z), 0)
    elseif ngens == nothing
        return _genmat_fallback_generic(Z, gens)
    else
        return _genmat_fallback_ngens(Z, gens, ngens)
    end
end

function _genmat_fallback_generic(Z::AbstractZonotope{N}, gens) where {N}
    Gv = Vector{Vector{N}}()
    @inbounds for (i, g) in enumerate(gens)
        push!(Gv, g)
    end
    G = Matrix{N}(undef, dim(Z), length(Gv))
    @inbounds for (i, g) in enumerate(Gv)
        G[:, i] = g
    end
    return G
end

function _genmat_fallback_ngens(Z::AbstractZonotope{N}, gens, ngens) where {N}
    G = Matrix{N}(undef, dim(Z), ngens)
    @inbounds for (i, g) in enumerate(gens)
        G[:, i] = g
    end
    return G
end

# iterator that wraps the generator matrix
struct FallbackGeneratorIterator{M<:AbstractMatrix}
    G::M
    n_plus_one::Int

    FallbackGeneratorIterator(G::M) where {M<:AbstractMatrix} =
        new{M}(G, size(G, 2) + 1)
end

Base.length(it::FallbackGeneratorIterator) = it.n_plus_one - 1

Base.eltype(::Type{<:FallbackGeneratorIterator{<:AbstractMatrix{N}}}) where {N} =
    AbstractVector{N}

Base.eltype(::Type{<:FallbackGeneratorIterator{<:Matrix{N}}}) where {N} =
    Vector{N}

Base.eltype(::Type{<:FallbackGeneratorIterator{<:SparseMatrixCSC{N}}}) where {N} =
    SparseVector{N}

function Base.iterate(it::FallbackGeneratorIterator, state::Int=1)
    if state == it.n_plus_one
        return nothing
    end
    g = it.G[:, state]
    state += 1
    return (g, state)
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
    return FallbackGeneratorIterator(genmat(Z))
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


# --- common AbstractZonotope functions ---


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

Return the order of a zonotope.

### Input

- `Z` -- zonotope

### Output

A rational number representing the order of the zonotope.

### Notes

The order of a zonotope is defined as the quotient of its number of generators
and its dimension.
"""
function order(Z::AbstractZonotope)
    return ngens(Z) // dim(Z)
end


# --- ConvexSet interface functions ---


"""
    ρ(d::AbstractVector, Z::AbstractZonotope)

Return the support function of a zonotopic set in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotopic set

### Output

The support function of the zonotopic set in the given direction.

### Algorithm

The support value is ``cᵀ d + ‖Gᵀ d‖₁`` where ``c`` is the center and ``G`` is
the generator matrix of `Z`.
"""
function ρ(d::AbstractVector, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    return dot(c, d) + abs_sum(d, G)
end

"""
    σ(d::AbstractVector, Z::AbstractZonotope)

Return the support vector of a zonotopic set in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotopic set

### Output

A support vector in the given direction.
If the direction has norm zero, the vertex with ``ξ_i = 1 \\ \\ ∀ i = 1,…, p``
is returned.
"""
function σ(d::AbstractVector, Z::AbstractZonotope)
    G = genmat(Z)
    return center(Z) .+ G * sign_cadlag.(At_mul_B(G, d))
end

"""
    ∈(x::AbstractVector, Z::AbstractZonotope; solver=nothing)

Check whether a given point is contained in a zonotopic set.

### Input

- `x`      -- point/vector
- `Z`      -- zonotopic set
- `solver` -- (optional, default: `nothing`) the backend used to solve the
              linear program

### Output

`true` iff ``x ∈ Z``.

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

The membership problem is computed by stating and solving the following linear
program.
Let ``p`` and ``n`` be the number of generators and ambient dimension,
respectively.
We consider the minimization of ``x_0`` in the ``p+1``-dimensional space of
elements ``(x_0, ξ_1, …, ξ_p)`` constrained to ``0 ≤ x_0 ≤ ∞``,
``ξ_i ∈ [-1, 1]`` for all ``i = 1, …, p``, and such that ``x-c = Gξ`` holds.
If a feasible solution exists, the optimal value ``x_0 = 0`` is achieved.
"""
function ∈(x::AbstractVector, Z::AbstractZonotope; solver=nothing)
    @assert length(x) == dim(Z)

    p, n = ngens(Z), dim(Z)
    N = promote_type(eltype(x), eltype(Z))
    # (n+1) x (p+1) matrix with block-diagonal blocks 1 and genmat(Z)
    A = [[one(N); zeros(N, p)]'; [zeros(N, n) genmat(Z)]]
    b = [zero(N); (x - center(Z))]
    lbounds = [zero(N); fill(-one(N), p)]
    ubounds = [N(Inf); ones(N, p)]
    sense = ['>'; fill('=', n)]
    obj = [one(N); zeros(N, p)]

    if solver == nothing
        solver = default_lp_solver(N)
    end
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    return is_lp_optimal(lp.status) # Infeasible or Unbounded => false
end

"""
    linear_map(M::AbstractMatrix, Z::AbstractZonotope)

Concrete linear map of a zonotopic set.

### Input

- `M` -- matrix
- `Z` -- zonotopic set

### Output

The zonotope obtained by applying the linear map to the center and generators
of ``Z``.
"""
function linear_map(M::AbstractMatrix, Z::AbstractZonotope)
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(Z))"

    c = M * center(Z)
    gi = M * genmat(Z)
    return Zonotope(c, gi)
end

"""
    translate(Z::AbstractZonotope, v::AbstractVector)

Translate (i.e., shift) a zonotope by a given vector.

### Input

- `Z`     -- zonotope
- `v`     -- translation vector

### Output

A translated zonotope.

### Notes

See also [`translate!(Z::AbstractZonotope, v::AbstractVector)`](@ref) for the
in-place version.

### Algorithm

We add the translation vector to the center of the zonotope.
"""
function translate(Z::AbstractZonotope, v::AbstractVector)
    return translate!(copy(Z), v)
end

"""
    translate!(Z::AbstractZonotope, v::AbstractVector)

Translate (i.e., shift) a zonotope by a given vector in-place.

### Input

- `Z`     -- zonotope
- `v`     -- translation vector

### Output

A translated zonotope.

### Notes

See also [`translate(Z::AbstractZonotope, v::AbstractVector)`](@ref) for the
out-of-place version.

### Algorithm

We add the translation vector to the center of the zonotope.
"""
function translate!(Z::AbstractZonotope, v::AbstractVector)
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(Z)
    c .+= v
    G = genmat(Z)
    return Z
end

# --- AbstractPolytope interface functions ---


"""
    vertices_list(Z::AbstractZonotope; [apply_convex_hull]::Bool=true)

Return the vertices of a zonotopic set.

### Input

- `Z`                 -- zonotopic set
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         computation with the convex hull of the points

### Output

List of vertices as a vector of vectors.

### Algorithm

#### Two-dimensional case

We use a trick to speed up enumerating vertices of 2-dimensional zonotopic
sets with all generators in the first quadrant or third quadrant (same sign).
Namely, sort the generators in angle and add them clockwise in increasing
order and anticlockwise in decreasing order, the algorithm detail:
https://math.stackexchange.com/q/3356460

To avoid cumulative sum from both directions separately, we build a 2d index matrix
to sum generators for both directions in one matrix-vector product.

#### General case

If the zonotopic set has ``p`` generators, each vertex is the result of summing
the center with some linear combination of generators, where the combination
factors are ``ξ_i ∈ \\{-1, 1\\}``.

There are at most ``2^p`` distinct vertices. Use the flag `apply_convex_hull` to
control whether a convex hull algorithm is applied to the vertices computed by
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
            return _vertices_list_2D_order_one_half(c, G, apply_convex_hull=apply_convex_hull)
        elseif p == 2
            return _vertices_list_2D_order_one(c, G, apply_convex_hull=apply_convex_hull)
        else
            return _vertices_list_2D(c, G, apply_convex_hull=apply_convex_hull)
        end

    else
        Gred = remove_zero_columns(G)
        return _vertices_list_iterative(c, Gred, apply_convex_hull=apply_convex_hull)
    end
end

"""
    constraints_list(P::AbstractZonotope)

Return the list of constraints defining a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

The list of constraints of the zonotopic set.

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
    constraints_list(Z::AbstractZonotope{N}) where {N<:AbstractFloat}

Return the list of constraints defining a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

The list of constraints of the zonotopic set.

### Notes

The algorithm assumes that no generator is redundant.
The result has ``2 \\binom{p}{n-1}`` (with ``p`` being the number of generators
and ``n`` being the ambient dimension) constraints, which is optimal under this
assumption.

If ``p < n`` or the generator matrix is not full rank, we fall back to the
(slower) computation based on the vertex representation.

### Algorithm

We follow the algorithm presented in *Althoff, Stursberg, Buss: Computing
Reachable Sets of Hybrid Systems Using a Combination of Zonotopes and Polytopes.
2009.*

The one-dimensional case is not covered by that algorithm; we manually handle
this case.
"""
function constraints_list(Z::AbstractZonotope{N}) where {N<:AbstractFloat}
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

    G = genmat(Z)
    c = center(Z)
    m = binomial(p, n - 1)
    constraints = Vector{LinearConstraint{N, Vector{N}}}()
    sizehint!(constraints, 2m)
    for columns in StrictlyIncreasingIndices(p, n-1)
        c⁺ = cross_product(view(G, :, columns))
        iszero(c⁺) && continue
        normalize!(c⁺, 2)

        Δd = sum(abs, transpose(G) * c⁺)

        d⁺ = dot(c⁺, c) + Δd
        c⁻ = -c⁺
        d⁻ = -d⁺ + 2 * Δd  # identical to dot(c⁻, c) + Δd

        push!(constraints, LinearConstraint(c⁺, d⁺))
        push!(constraints, LinearConstraint(c⁻, d⁻))
    end
    return constraints
end

function _constraints_list_1d(Z::AbstractZonotope{N}) where {N}
    c = center(Z, 1)
    g = sum(abs, genmat(Z))
    return [LinearConstraint([N(1)], c + g), LinearConstraint([N(-1)], g - c)]
end

"""
    split(Z::AbstractZonotope, j::Int)

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
function split(Z::AbstractZonotope, j::Int)
    return _split(convert(Zonotope, Z), j)
end

"""
    split(Z::AbstractZonotope, gens::AbstractVector{Int}, nparts::AbstractVector{Int})

Split a zonotope along the given generators into a vector of zonotopes.

### Input

- `Z`    -- zonotope
- `gens` -- vector of indices of the generators to be split
- `n`    -- vector of integers describing the number of partitions in the
            corresponding generator

### Output

The zonotopes obtained by splitting `Z` into `2^{n_i}` zonotopes for each
generator `i` such that their union is `Z` and their intersection is
possibly non-empty.

### Examples

Splitting of a two-dimensional zonotope along its first generator:

```jldoctest zonotope_label
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> split(Z, [1], [1])
2-element Vector{Zonotope{Float64, Vector{Float64}, Matrix{Float64}}}:
 Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([0.95, 0.0], [0.05 0.0; 0.0 0.1])
 Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.05, 0.0], [0.05 0.0; 0.0 0.1])
```
Here, the first vector in the arguments corresponds to the zonotope's
generator to be split, and the second vector corresponds to the exponent of
`2^n` parts that the zonotope will be split into along the corresponding generator.

Splitting of a two-dimensional zonotope along its generators:

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
Here the zonotope is split along both of its generators, each time into four parts.
"""
function split(Z::AbstractZonotope, gens::AbstractVector{Int}, nparts::AbstractVector{Int})
    return _split(convert(Zonotope, Z), gens, nparts)
end

function project(Z::AbstractZonotope{N}, block::AbstractVector{Int}; kwargs...) where {N}
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

Remove all redundant (pairwise linearly dependent) generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

A new zonotope with fewer generators, or the same zonotope if no generator could
be removed.

### Algorithm

By default this function returns the input zonotope. Subtypes of
`AbstractZonotope` where generators can be removed have to define a new method.
"""
function remove_redundant_generators(Z::AbstractZonotope)
    return Z  # fallback implementation
end

# ==================================
# Zonotope order reduction methods
# ==================================

"""
    AbstractReductionMethod

Abstract supertype for zonotope order reduction methods.
"""
abstract type AbstractReductionMethod end

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

"""
    reduce_order(Z::AbstractZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())

Reduce the order of a zonotope by overapproximating with a zonotope with fewer
generators.

### Input

- `Z`      -- zonotope
- `r`      -- desired order
- `method` -- (optional, default: `GIR05()`) the reduction method used

### Output

A new zonotope with fewer generators, if possible.

### Algorithm

The available algorithms are:

```jldoctest; setup = :(using LazySets: subtypes, AbstractReductionMethod)
julia> subtypes(AbstractReductionMethod)
2-element Vector{Any}:
 LazySets.COMB03
 LazySets.GIR05
```

See the documentation of each algorithm for references. These methods split the
given zonotope `Z` into two zonotopes, `K` and `L`, where `K` contains the most
"representative" generators and `L` contains the generators that are reduced,
`Lred`, using a box overapproximation. We follow the notation from [YS18]. See
also [KSA17].

- [KSA17] Kopetzki, A. K., Schürmann, B., & Althoff, M. (2017). *Methods for
order reduction of zonotopes.* CDC.
- [YS18] Yang, X., & Scott, J. K. (2018). *A comparison of zonotope order
reduction techniques.* Automatica, 95.
"""
function reduce_order(Z::AbstractZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    r >= 1 || throw(ArgumentError("the target order should be at least 1, " *
        "but it is $r"))
    n = dim(Z)
    p = ngens(Z)

    # if r is bigger than the order of Z => don't reduce
    (r * n >= p) && return Z

    c = center(Z)
    G = genmat(Z)

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
