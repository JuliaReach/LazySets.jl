import Base: ∈, split

export AbstractZonotope,
       genmat,
       generators,
       ngens,
       order,
       togrep,
       split!

"""
    AbstractZonotope{N<:Real} <: AbstractCentrallySymmetricPolytope{N}

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
4-element Array{Any,1}:
 AbstractHyperrectangle
 HParallelotope
 LineSegment
 Zonotope
```
"""
abstract type AbstractZonotope{N<:Real} <: AbstractCentrallySymmetricPolytope{N} end

isconvextype(::Type{<:AbstractZonotope}) = true

# --- fallback implementations of AbstractZonotope functions ---


"""
    genmat_fallback(Z::AbstractZonotope{N}) where {N<:Real}

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
                         ngens=nothing) where {N<:Real}
    if isempty(gens)
        return Matrix{N}(undef, dim(Z), 0)
    elseif ngens == nothing
        return _genmat_fallback_generic(Z, gens)
    else
        return _genmat_fallback_ngens(Z, gens, ngens)
    end
end

function _genmat_fallback_generic(Z::AbstractZonotope{N}, gens) where {N<:Real}
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

function _genmat_fallback_ngens(Z::AbstractZonotope{N}, gens, ngens) where {N<:Real}
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
    generators_fallback(Z::AbstractZonotope{N}) where {N<:Real}

Fallback definition of `generators` for zonotopic sets.

### Input

- `Z` -- zonotopic set

### Output

An iterator over the generators of `Z`.
"""
function generators_fallback(Z::AbstractZonotope{N}) where {N<:Real}
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


# --- LazySet interface functions ---


"""
    ρ(d::AbstractVector{N}, Z::AbstractZonotope{N}) where {N<:Real}

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
function ρ(d::AbstractVector{N}, Z::AbstractZonotope{N}) where {N<:Real}
    c = center(Z)
    G = genmat(Z)
    return dot(c, d) + _abs_sum(d, G)
end

"""
    σ(d::AbstractVector{N}, Z::AbstractZonotope{N}) where {N<:Real}

Return the support vector of a zonotopic set in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotopic set

### Output

A support vector in the given direction.
If the direction has norm zero, the vertex with ``ξ_i = 1 \\ \\ ∀ i = 1,…, p``
is returned.
"""
function σ(d::AbstractVector{N}, Z::AbstractZonotope{N}) where {N<:Real}
    G = genmat(Z)
    return center(Z) .+ G * sign_cadlag.(_At_mul_B(G, d))
end

"""
    ∈(x::AbstractVector{N}, Z::AbstractZonotope{N};
      solver=default_lp_solver(N)) where {N<:Real}

Check whether a given point is contained in a zonotopic set.

### Input

- `x`      -- point/vector
- `Z`      -- zonotopic set
- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

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
function ∈(x::AbstractVector{N}, Z::AbstractZonotope{N};
           solver=default_lp_solver(N)) where {N<:Real}
    @assert length(x) == dim(Z)

    p, n = ngens(Z), dim(Z)
    # (n+1) x (p+1) matrix with block-diagonal blocks 1 and genmat(Z)
    A = [[one(N); zeros(N, p)]'; [zeros(N, n) genmat(Z)]]
    b = [zero(N); (x - center(Z))]
    lbounds = [zero(N); fill(-one(N), p)]
    ubounds = [N(Inf); ones(N, p)]
    sense = ['>'; fill('=', n)]
    obj = [one(N); zeros(N, p)]

    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    return (lp.status == :Optimal) # Infeasible or Unbounded => false
end

"""
    linear_map(M::AbstractMatrix{N}, Z::AbstractZonotope{N}) where {N<:Real}

Concrete linear map of a zonotopic set.

### Input

- `M` -- matrix
- `Z` -- zonotopic set

### Output

The zonotope obtained by applying the linear map to the center and generators
of ``Z``.
"""
function linear_map(M::AbstractMatrix{N}, Z::AbstractZonotope{N}
                   ) where {N<:Real}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(Z))"

    c = M * center(Z)
    gi = M * genmat(Z)
    return Zonotope(c, gi)
end

"""
    translate(Z::AbstractZonotope{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a zonotope by a given vector.

### Input

- `Z`     -- zonotope
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated zonotope.

### Notes

The generator matrix is shared with the original zonotope if `share == true`.

### Algorithm

We add the vector to the center of the zonotope.
"""
function translate(Z::AbstractZonotope{N}, v::AbstractVector{N};
                   share::Bool=false) where {N<:Real}
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(Z) + v
    G = share ? genmat(Z) : copy(genmat(Z))
    return Zonotope(c, G)
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(Z::AbstractZonotope{N}; [apply_convex_hull]::Bool=true
                 ) where {N<:Real}

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
function vertices_list(Z::AbstractZonotope{N};
                       apply_convex_hull::Bool=true) where {N<:Real}
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
    constraints_list(P::AbstractZonotope{N}) where {N<:Real}

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
function constraints_list(Z::AbstractZonotope{N}) where {N<:Real}
    return constraints_list(VPolytope(vertices_list(Z)))
end

"""
    constraints_list(Z::AbstractZonotope{N}; check_full_rank::Bool=true
                    ) where {N<:AbstractFloat}

Return the list of constraints defining a zonotopic set.

### Input

- `Z`               -- zonotopic set
- `check_full_rank` -- (optional; default: `true`) flag for checking whether the
                       generator matrix has full rank

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
this case, assuming that there is only one generator.
"""
function constraints_list(Z::AbstractZonotope{N}; check_full_rank::Bool=true
                         ) where {N<:AbstractFloat}
    G = genmat(Z)
    p = ngens(Z)
    n = dim(Z)

    # use fallback implementation if order < 1 or matrix is not full rank
    if p < n || (check_full_rank && rank(G) < n)
        return invoke(constraints_list, Tuple{AbstractZonotope{<:Real}}, Z)
    end

    # special handling of 1D case
    if n == 1
        c = center(Z, 1)
        g = sum(abs.(view(G, 1, :)))
        constraints = [LinearConstraint([N(1)], c + g),
                       LinearConstraint([N(-1)], g - c)]
        return constraints
    end

    c = center(Z)
    m = binomial(p, n - 1)
    constraints = Vector{LinearConstraint{N, Vector{N}}}()
    sizehint!(constraints, 2m)
    for columns in StrictlyIncreasingIndices(p, n-1)
        c⁺ = cross_product(view(G, :, columns))
        iszero(c⁺) && continue
        normalize!(c⁺, 2)

        Δd = sum(abs.(transpose(G) * c⁺))

        d⁺ = dot(c⁺, c) + Δd
        c⁻ = -c⁺
        d⁻ = -d⁺ + 2 * Δd  # identical to dot(c⁻, c) + Δd

        push!(constraints, LinearConstraint(c⁺, d⁺))
        push!(constraints, LinearConstraint(c⁻, d⁻))
    end
    return constraints
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
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> split(Z, [1], [1])
2-element Array{Zonotope{Float64,Array{Float64,1},Array{Float64,2}},1}:
 Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.95, 0.0], [0.05 0.0; 0.0 0.1])
 Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.05, 0.0], [0.05 0.0; 0.0 0.1])
```
Here, the first vector in the arguments corresponds to the zonotope's
generator to be split, and the second vector corresponds to the exponent of
`2^n` parts that the zonotope will be split into along the corresponding generator.

Splitting of a two-dimensional zonotope along its generators:

```
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> split(Z, [1, 2], [2, 2])
16-element Array{Zonotope{Float64,Array{Float64,1},Array{Float64,2}},1}:
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.925, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.925, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.925, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.925, 0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.975, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.975, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.975, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([0.975, 0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.025, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.025, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.025, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.025, 0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.075, -0.075], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.075, -0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.075, 0.025], [0.025 0.0; 0.0 0.025])
Zonotope{Float64,Array{Float64,1},Array{Float64,2}}([1.075, 0.075], [0.025 0.0; 0.0 0.025])
```
Here the zonotope is split along both of its generators, each time into four parts.
"""
function split(Z::AbstractZonotope, gens::AbstractVector{Int}, nparts::AbstractVector{Int})
    return _split(convert(Zonotope, Z), gens, nparts)
end
