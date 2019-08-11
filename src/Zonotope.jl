import Base: rand,
             ∈,
             split

export Zonotope,
       order,
       linear_map,
       scale,
       ngens,
       reduce_order,
       constraints_list,
       generators

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
        if remove_zero_generators
            generators = delete_zero_columns!(generators)
        end
        new{N}(center, generators)
    end
end

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


# --- AbstractPolytope interface functions ---


"""
    vertices_list(Z::Zonotope{N};
                  [apply_convex_hull]::Bool=true)::Vector{Vector{N}} where {N<:Real}

Return the vertices of a zonotope.

### Input

- `Z`                 -- zonotope
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         computation with the convex hull of the points

### Output

List of vertices as a vector of vectors.

### Algorithm

If the zonotope has ``p`` generators, each vertex is the result of summing the
center with some linear combination of generators, where the combination
factors are ``ξ_i ∈ \\{-1, 1\\}``.

There are at most ``2^p`` distinct vertices. Use the flag `apply_convex_hull` to
control whether a convex hull algorithm is applied to the vertices computed by this
method; otherwise, redundant vertices may be present.
"""
function vertices_list(Z::Zonotope{N};
                       apply_convex_hull::Bool=true)::Vector{Vector{N}} where {N<:Real}
    p = ngens(Z)
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, 2^p)

    for ξi in Iterators.product([[1, -1] for i = 1:p]...)
        push!(vlist, Z.center .+ Z.generators * collect(ξi))
    end

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

# --- LazySet interface functions ---

"""
    ρ(d::AbstractVector{N}, Z::Zonotope{N}) where {N<:Real}

Return the support function of a zonotope in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotope

### Output

The support function of the zonotope in the given direction.

### Algorithm

The support value is ``cᵀ d + ‖Gᵀ d‖₁`` where ``c`` is the center and ``G`` is
the generator matrix of `Z`.

"""
function ρ(d::AbstractVector{N}, Z::Zonotope{N}) where {N<:Real}
    return dot(center(Z), d) + sum(abs.(transpose(Z.generators) * d))
end

"""
    σ(d::AbstractVector{N}, Z::Zonotope{N}) where {N<:Real}

Return the support vector of a zonotope in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotope

### Output

Support vector in the given direction.
If the direction has norm zero, the vertex with ``ξ_i = 1 \\ \\ ∀ i = 1,…, p``
is returned.
"""
function σ(d::AbstractVector{N}, Z::Zonotope{N}) where {N<:Real}
    return Z.center .+ Z.generators * sign_cadlag.(_At_mul_B(Z.generators, d))
end

"""
    ∈(x::AbstractVector{N}, Z::Zonotope{N};
      solver=GLPKSolverLP(method=:Simplex))::Bool where {N<:Real}

Check whether a given point is contained in a zonotope.

### Input

- `x`      -- point/vector
- `Z`      -- zonotope
- `solver` -- (optional, default: `GLPKSolverLP(method=:Simplex)`) the backend
              used to solve the linear program

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
program with the simplex method.
Let ``p`` and ``n`` be the number of generators and ambient dimension,
respectively.
We consider the minimization of ``x_0`` in the ``p+1``-dimensional space of
elements ``(x_0, ξ_1, …, ξ_p)`` constrained to ``0 ≤ x_0 ≤ ∞``,
``ξ_i ∈ [-1, 1]`` for all ``i = 1, …, p``, and such that ``x-c = Gξ`` holds.
If a feasible solution exists, the optimal value ``x_0 = 0`` is achieved.

### Notes

This function is parametric in the number type `N`. For exact arithmetic use
an appropriate backend, e.g. `solver=GLPKSolverLP(method=:Exact)`.
"""
function ∈(x::AbstractVector{N}, Z::Zonotope{N};
           solver=GLPKSolverLP(method=:Simplex))::Bool where {N<:Real}
    @assert length(x) == dim(Z)

    p, n = ngens(Z), dim(Z)
    # (n+1) x (p+1) matrix with block-diagonal blocks 1 and Z.generators
    A = [[one(N); zeros(N, p)]'; [zeros(N, n) Z.generators]]
    b = [zero(N); (x - Z.center)]
    lbounds = [zero(N); fill(-one(N), p)]
    ubounds = [N(Inf); ones(N, p)]
    sense = ['>'; fill('=', n)]
    obj = [one(N); zeros(N, p)]

    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)
    return (lp.status == :Optimal) # Infeasible or Unbounded => false
end

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
    order(Z::Zonotope)::Rational

Return the order of a zonotope.

### Input

- `Z` -- zonotope

### Output

A rational number representing the order of the zonotope.

### Notes

The order of a zonotope is defined as the quotient of its number of generators
and its dimension.
"""
function order(Z::Zonotope)::Rational
    return ngens(Z) // dim(Z)
end

"""
    linear_map(M::AbstractMatrix{N}, Z::Zonotope{N}) where {N<:Real}

Concrete linear map of a zonotope.

### Input

- `M` -- matrix
- `Z` -- zonotope

### Output

The zonotope obtained by applying the linear map to the center and generators
of ``Z``.
"""
function linear_map(M::AbstractMatrix{N}, Z::Zonotope{N}) where {N<:Real}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(Z))"

    c = M * Z.center
    gi = M * Z.generators
    return Zonotope(c, gi)
end

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
    ngens(Z::Zonotope)::Int

Return the number of generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

Integer representing the number of generators.
"""
ngens(Z::Zonotope)::Int = size(Z.generators, 2)

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

"""
    constraints_list(P::Zonotope{N}) where {N<:Real}

Return the list of constraints defining a zonotope.

### Input

- `Z` -- zonotope

### Output

The list of constraints of the zonotope.

### Algorithm

This is the (inefficient) fallback implementation for rational numbers.
It first computes the vertices and then converts the corresponding polytope
to constraint representation.
"""
function constraints_list(Z::Zonotope{N}) where {N<:Real}
    return constraints_list(VPolytope(vertices_list(Z)))
end

"""
    constraints_list(Z::Zonotope{N}) where {N<:AbstractFloat}

Return the list of constraints defining a zonotope.

### Input

- `Z`               -- zonotope
- `check_full_rank` -- (optional; default: `true`) flag for checking whether the
                       generator matrix has full rank

### Output

The list of constraints of the zonotope.

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
function constraints_list(Z::Zonotope{N}; check_full_rank::Bool=true
                         ) where {N<:AbstractFloat}
    G = genmat(Z)
    p = ngens(Z)
    n = dim(Z)

    # use fallback implementation if order < 1 or matrix is not full rank
    if p < n || (check_full_rank && rank(G) < n)
        return invoke(constraints_list, Tuple{Zonotope{<:Real}}, Z)
    end

    # special handling of 1D case
    if n == 1
        if p > 1
            error("1D-zonotope constraints currently only support a single " *
                  "generator")
        end

        c = Z.center[1]
        g = G[:, 1][1]
        constraints = [LinearConstraint([N(1)], c + g),
                       LinearConstraint([N(-1)], g - c)]
        return constraints
    end

    i = 0
    c = Z.center
    m = binomial(p, n - 1)
    constraints = Vector{LinearConstraint{N, Vector{N}}}(undef, 2 * m)
    for columns in StrictlyIncreasingIndices(p, n-1)
        i += 1
        c⁺ = cross_product(view(G, :, columns))
        normalize!(c⁺, 2)

        Δd = sum(abs.(transpose(G) * c⁺))

        d⁺ = dot(c⁺, c) + Δd
        c⁻ = -c⁺
        d⁻ = -d⁺ + 2 * Δd  # identical to dot(c⁻, c) + Δd

        constraints[i] = LinearConstraint(c⁺, d⁺)
        constraints[i + m] = LinearConstraint(c⁻, d⁻)
    end
    @assert i == m "expected 2*$m constraints, but only created 2*$i"
    return constraints
end

"""
    translate(Z::Zonotope{N}, v::AbstractVector{N}; share::Bool=false
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
function translate(Z::Zonotope{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(Z) + v
    generators = share ? Z.generators : copy(Z.generators)
    return Zonotope(c, generators)
end
