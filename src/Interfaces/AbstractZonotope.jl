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

- `genmat(::AbstractZonotope)` -- return the generator matrix

- `generators(::AbstractZonotope)` -- return an iterator over the generators

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
 RotatedHyperrectangle
 Zonotope
```
"""
abstract type AbstractZonotope{N} <: AbstractCentrallySymmetricPolytope{N} end

"""
    genmat_fallback(Z::AbstractZonotope{N};
                    [gens]=generators(Z), [ngens]=nothing) where {N}

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
function genmat_fallback(Z::AbstractZonotope{N};
                         gens=generators(Z), ngens=nothing) where {N}
    if isempty(gens)
        return Matrix{N}(undef, dim(Z), 0)
    elseif isnothing(ngens)
        return _genmat_fallback_generic(Z, gens)
    else
        return _genmat_fallback_ngens(Z, gens, ngens)
    end
end

function _genmat_fallback_generic(Z::AbstractZonotope{N}, gens) where {N}
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
    ρ(d::AbstractVector, Z::AbstractZonotope)

Evaluate the support function of a zonotopic set in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotopic set

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support value is ``cᵀ d + ‖Gᵀ d‖₁``, where ``c`` is the center and ``G`` is
the generator matrix of `Z`.
"""
function ρ(d::AbstractVector, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    return dot(c, d) + abs_sum(d, G)
end

"""
    σ(d::AbstractVector, Z::AbstractZonotope)

Return a support vector of a zonotopic set in a given direction.

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

The membership problem is reduced to the following linear program.
Let ``p`` and ``n`` be the number of generators and ambient dimension,
respectively.
We consider the ``p``-dimensional space of elements ``(ξ_1, …, ξ_p)``
constrained to ``ξ_i ∈ [-1, 1]`` for all ``i = 1, …, p`` such that
``x-c = Gξ`` holds.
"""
function ∈(x::AbstractVector, Z::AbstractZonotope; solver=nothing)
    @assert length(x) == dim(Z)

    p, n = ngens(Z), dim(Z)
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
    return is_lp_optimal(lp.status) # Infeasible or Unbounded => false
end

"""
    linear_map(M::AbstractMatrix, Z::AbstractZonotope)

Concrete linear map of a zonotopic set.

### Input

- `M` -- matrix
- `Z` -- zonotopic set

### Output

The zonotope obtained by applying the linear map.

### Algorithm

We apply the linear map to the center and the generators.
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

Translate (i.e., shift) a zonotopic set by a given vector.

### Input

- `Z` -- zonotopic set
- `v` -- translation vector

### Output

A translated zonotopic set.

### Notes

See also [`translate!(Z::AbstractZonotope, v::AbstractVector)`](@ref) for the
in-place version.
"""
function translate(Z::AbstractZonotope, v::AbstractVector)
    return translate!(copy(Z), v)
end

"""
    translate!(Z::AbstractZonotope, v::AbstractVector)

Translate (i.e., shift) a zonotopic set by a given vector in-place.

### Input

- `Z` -- zonotopic set
- `v` -- translation vector

### Output

A translated zonotopic set.

### Notes

See also [`translate(Z::AbstractZonotope, v::AbstractVector)`](@ref) for the
out-of-place version.

### Algorithm

We add the translation vector to the center of the zonotopic set.
"""
function translate!(Z::AbstractZonotope, v::AbstractVector)
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(Z)
    c .+= v
    return Z
end

"""
    vertices_list(Z::AbstractZonotope; [apply_convex_hull]::Bool=true)

Return a list of the vertices of a zonotopic set.

### Input

- `Z`                 -- zonotopic set
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         computation with the convex hull of the points

### Output

A list of the vertices.

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
        return _vertices_list_zonotope_2D_positive(c, G; apply_convex_hull=apply_convex_hull)
    else
        # FIXME generalized 2D vertices list function is not implemented yet
        # See LazySets#2209
        return _vertices_list_zonotope_iterative(c, G; apply_convex_hull=apply_convex_hull)
    end
end

function _vertices_list_zonotope_2D_positive(c::AbstractVector{N}, G::AbstractMatrix{N};
                                             apply_convex_hull::Bool) where {N}
    n, p = size(G)

    # TODO special case p = 1 or p = 2 ?

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

function _vertices_list_zonotope_iterative(c::VN, G::MN;
                                           apply_convex_hull::Bool) where {N,VN<:AbstractVector{N},
                                                                           MN<:AbstractMatrix{N}}
    p = size(G, 2)
    vlist = Vector{VN}()
    sizehint!(vlist, 2^p)

    for ξi in Iterators.product([(1, -1) for i in 1:p]...)
        push!(vlist, c .+ G * collect(ξi))
    end

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

# special case 2D zonotope of order 1/2
function _vertices_list_zonotope_2D_order_one_half(c::VN, G::MN;
                                                   apply_convex_hull::Bool) where {N,
                                                                                   VN<:AbstractVector{N},
                                                                                   MN}
    vlist = Vector{VN}(undef, 2)
    g = view(G, :, 1)
    @inbounds begin
        vlist[1] = c .+ g
        vlist[2] = c .- g
    end
    return apply_convex_hull ? _two_points_2d!(vlist) : vlist
end

# special case 2D zonotope of order 1
function _vertices_list_zonotope_2D_order_one(c::VN, G::MN;
                                              apply_convex_hull::Bool) where {N,
                                                                              VN<:AbstractVector{N},
                                                                              MN}
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
    constraints_list(P::AbstractZonotope)

Return a list of constraints defining a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

A list of constraints of the zonotopic set.

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

Return a list of constraints defining a zonotopic set.

### Input

- `Z` -- zonotopic set

### Output

A list of constraints of the zonotopic set.

### Notes

The main algorithm assumes that the generator matrix is full rank.
The result has ``2 \\binom{p}{n-1}`` (with ``p`` being the number of generators
and ``n`` being the ambient dimension) constraints, which is optimal under this
assumption.
If this assumption is not given, the implementation tries to work around.

### Algorithm

We follow the algorithm presented in [1]. Three cases are not covered by that
algorithm, so we handle them separately. The first case is zonotopes in one
dimension. The second case is that there are fewer generators than dimensions,
``p < n``, or the generator matrix is not full rank, in which case we fall back
to the (slower) computation based on the vertex representation. The third case
is that the zonotope is flat in some dimensions, in which case we project the
zonotope to the non-flat dimensions and extend the result later.

[1] Althoff, Stursberg, Buss. *Computing Reachable Sets of Hybrid Systems Using
a Combination of Zonotopes and Polytopes*. 2009.
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

This function implements [Prop. 3, 1], which we state next. The zonotopic set
``Z = ⟨c, g^{(1, …, p)}⟩`` is split into:

```math
Z₁ = ⟨c - \\frac{1}{2}g^{(j)}, (g^{(1, …,j-1)}, \\frac{1}{2}g^{(j)}, g^{(j+1, …, p)})⟩ \\\\
Z₂ = ⟨c + \\frac{1}{2}g^{(j)}, (g^{(1, …,j-1)}, \\frac{1}{2}g^{(j)}, g^{(j+1, …, p)})⟩,
```
such that ``Z₁ ∪ Z₂ = Z`` and ``Z₁ ∩ Z₂ = Z^*``, where

```math
Z^* = ⟨c, (g^{(1,…,j-1)}, g^{(j+1,…, p)})⟩.
```

[1] Althoff, M., Stursberg, O., & Buss, M. *Reachability analysis of nonlinear
systems with uncertain parameters using conservative linearization*. CDC 2008.
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
function project(Z::AbstractZonotope{N}, block::AbstractVector{Int};
                 kwargs...) where {N}
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

"""
    AbstractReductionMethod

Abstract supertype for order-reduction methods of a zonotopic set.
"""
abstract type AbstractReductionMethod end

"""
    GIR05 <: AbstractReductionMethod

Zonotope order-reduction method from [1].

- [1] A. Girard. *Reachability of Uncertain Linear Systems Using Zonotopes*.
HSCC 2005.
"""
struct GIR05 <: AbstractReductionMethod end

"""
    COMB03 <: AbstractReductionMethod

Zonotope order-reduction method from [1].

- [1] C. Combastel. *A state bounding observer based on zonotopes*. ECC 2003.
"""
struct COMB03 <: AbstractReductionMethod end

"""
    ASB10 <: AbstractReductionMethod

Zonotope order-reduction method from [1].

- [1] Althoff, M., Stursberg, O., & Buss, M. *Computing reachable sets of hybrid
systems using a combination of zonotopes and polytopes*. Nonlinear analysis:
hybrid systems 2010.
"""
struct ASB10 <: AbstractReductionMethod end

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
3-element Vector{Any}:
 LazySets.ASB10
 LazySets.COMB03
 LazySets.GIR05
```

See the documentation of each algorithm for references. These methods split the
given zonotopic set `Z` into two zonotopes, `K` and `L`, where `K` contains the
most "representative" generators and `L` contains the generators that are
reduced, `Lred`, using a box overapproximation. We follow the notation from [1].
See also [2].

- [1] Yang, X., & Scott, J. K. *A comparison of zonotope order reduction
techniques*. Automatica 2018.
- [2] Kopetzki, A. K., Schürmann, B., & Althoff, M. *Methods for order reduction
of zonotopes*. CDC 2017.
- [3] Althoff, M., Stursberg, O., & Buss, M. *Computing reachable sets of hybrid
systems using a combination of zonotopes and polytopes*. Nonlinear analysis:
hybrid systems 2010.
"""
function reduce_order(Z::AbstractZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    r >= 1 || throw(ArgumentError("the target order should be at least 1, " *
                                  "but it is $r"))
    n = dim(Z)
    p = ngens(Z)

    # if r is bigger than the order of Z => do not reduce
    (r * n >= p) && return Z

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

    # compute interval hull of L
    Lred = _approximate_reduce_order(c, G, view(indices, (m + 1):p), method)

    # concatenate non-reduced and reduced generators
    Gred = _hcat_KLred(G, view(indices, 1:m), Lred)

    return Zonotope(c, Gred)
end

# approximate with a box
function _approximate_reduce_order(c, G, indices, method::Union{COMB03,GIR05})
    return _interval_hull(G, indices)
end

# approximate with a parallelotope
function _approximate_reduce_order(c, G, indices, method::ASB10)
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
# ||⋅||₁ - ||⋅||∞ difference. The generator index with highest score goes first.
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
    reflect(Z::AbstractZonotope)

Concrete reflection of a zonotopic set `Z`, resulting in the reflected set `-Z`.

### Input

- `Z` -- zonotopic set

### Output

A `Zonotope` representing `-Z`.

### Algorithm

If ``Z`` has center ``c`` and generator matrix ``G``, then ``-Z`` has center
``-c`` and generator matrix ``G``. For the latter, observe that ``G`` and ``-G``
behave the same way.
"""
function reflect(Z::AbstractZonotope)
    return Zonotope(-center(Z), genmat(Z))
end
