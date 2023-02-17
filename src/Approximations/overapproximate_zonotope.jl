"""
    overapproximate(Z::Zonotope, ::Type{<:Zonotope}, r::Union{Integer, Rational})

Reduce the order of a zonotope.

### Input

- `Z`        -- zonotope
- `Zonotope` -- target set type
- `r`        -- desired order

### Output

A new zonotope with `r` generators, if possible.

### Algorithm

This method falls back to `reduce_order` with the default algorithm.
"""
function overapproximate(Z::Zonotope, ::Type{<:Zonotope},
                         r::Union{Integer, Rational})
    reduce_order(Z, r)
end

"""
    overapproximate(X::ConvexHull{N, <:AbstractZonotope, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the convex hull of two zonotopic sets.

### Input

- `X`         -- convex hull of two zonotopic sets
- `Zonotope`  -- target set type
- `algorithm` -- (optional; default: `"mean"`) choice of algorithm; possible
                 values are `"mean"` and `"join"`

### Output

A zonotope ``Z`` such that ``X ⊆ Z``.

### Algorithm

The algorithm can be controlled by the parameter `algorithm`.
Note that the results of the two implemented algorithms are generally
incomparable.

##### 'mean' method

If `algorithm == "mean"`, we choose the method proposed in [1].
The convex hull of two zonotopic sets ``Z₁`` and ``Z₂`` of the same order,
which we write

```math
Z_j = ⟨c^{(j)}, g^{(j)}_1, …, g^{(j)}_p⟩
```
for ``j = 1, 2``, can be overapproximated as follows:

```math
CH(Z_1, Z_2) ⊆ \\frac{1}{2}⟨c^{(1)}+c^{(2)}, g^{(1)}_1+g^{(2)}_1, …,
g^{(1)}_p+g^{(2)}_p, c^{(1)}-c^{(2)}, g^{(1)}_1-g^{(2)}_1, …, g^{(1)}_p-g^{(2)}_p⟩.
```

If the zonotope order is not the same, this algorithm calls
`reduce_order` to reduce the order to the minimum of the arguments.

It should be noted that the output zonotope is not necessarily the minimal
enclosing zonotope, which is in general expensive to obtain in high dimensions.
This is further investigated in [2].

##### 'join' method

If `algorithm == "join"`, we choose the method proposed in [3, Definition 1].
The convex hull ``X`` of two zonotopic sets ``Z₁`` and ``Z₂`` is
overapproximated by a zonotope ``Z₃`` such that the box approximation of ``X``
is identical with the box approximation of ``Z₃``.
Let ``□(X)`` denote the box approximation of ``X``.
The center of ``Z₃`` is the center of ``□(X)``.

The generator construction consists of two phases.
In the first phase, we construct generators ``g`` as a combination of one
generator from ``Z₁``, say, ``g₁``, with another generator from ``Z₂``, say,
``g₂``.
The entry of ``g`` in the ``i``-th dimension is given as

```math
    g[i] = \\arg\\min_{\\min(g₁[i], g₂[i]) ≤ x ≤ \\max(g₁[i], g₂[i])} |x|.
```

If ``g`` is the zero vector, it can be omitted.

In the second phase, we construct another generator for each dimension.
These generators are scaled unit vectors.
The following formula defines the sum of all those generators.

```math
    \\sup(□(X)) - c - ∑_g |g|
```

where ``c`` is the center of the new zonotope and the ``g``s are the generators
constructed in the first phase.

##### References

[1] *Reachability of Uncertain Linear Systems Using Zonotopes*. A. Girard.
    HSCC 2005.

[2] *Zonotopes as bounding volumes*. L. J. Guibas, A. T. Nguyen, L. Zhang.
    SODA 2003.

[3] *The zonotope abstract domain Taylor1+*. K. Ghorbal, E. Goubault, S. Putot.
    CAV 2009.
"""
function overapproximate(X::ConvexHull{N, <:AbstractZonotope, <:AbstractZonotope},
                         ::Type{<:Zonotope};
                         algorithm="mean") where {N}
    # execute specific algorithm
    if algorithm == "mean"
        return _overapproximate_convex_hull_zonotope_G05(X)
    elseif algorithm == "join"
        return _overapproximate_convex_hull_zonotope_GGP09(X)
    else
        error("algorithm $algorithm is not known")
    end
end

function _overapproximate_convex_hull_zonotope_G05(X::ConvexHull{N}) where {N}
    # reduce to the same order if possible
    m1, m2 = ngens(X.X), ngens(X.Y)
    if m1 < m2
        Z1 = X.X
        Z2 = reduce_order(X.Y, max(1, m1))
    elseif m1 > m2
        Z1 = reduce_order(X.X, max(1, m2))
        Z2 = X.Y
    else
        Z1 = X.X
        Z2 = X.Y
    end

    if order(Z2) > order(Z1)
        Z1, Z2 = Z2, Z1
    end

    c1 = center(Z1)
    c2 = center(Z2)
    G1 = genmat(Z1)
    G2 = genmat(Z2)
    c = (c1 + c2) / N(2)

    # the case of equal order is treated separately to avoid a slicing
    # (this creates a copy)
    if order(Z1) == order(Z2)
        G = hcat(G1 .+ G2,
                 c1 - c2,
                 G1 .- G2) / N(2)
    else
        G = hcat((G1[:, 1:ngens(Z2)] .+ G2) / N(2),
                 (c1 - c2) / N(2),
                 (G1[:, 1:ngens(Z2)] .- G2) / N(2),
                 G1[:, ngens(Z2)+1:end])
    end
    Z = Zonotope(c, G)
    return remove_zero_generators(Z)
end

function _overapproximate_convex_hull_zonotope_GGP09(X::ConvexHull{N}) where {N}
    Z1, Z2 = X.X, X.Y
    m = min(ngens(Z1), ngens(Z2))
    G1, G2 = genmat(Z1), genmat(Z2)
    n = dim(Z1)
    box = box_approximation(X)

    # new center: mid point of box approximation
    c = center(box)

    # the k-th new generator is a simple combination of the old k-th generators
    G = Vector{Vector{N}}()
    sizehint!(G, m + n)
    g_sum = zeros(N, n)
    @inbounds for j in 1:m
        g1, g2 = G1[:, j], G2[:, j]
        g = Vector{N}(undef, n)
        for i in 1:n
            gi_min, gi_max = g1[i] < g2[i] ? (g1[i], g2[i]) : (g2[i], g1[i])
            if gi_min <= zero(N) && gi_max >= zero(N)
                g[i] = zero(N)
            elseif abs(gi_min) <= abs(gi_max)
                g[i] = gi_min
            else
                g[i] = gi_max
            end
        end
        if !iszero(g)
            push!(G, g)
            g_sum += abs.(g)
        end
    end

    # one more new generator (a scaled unit vector) for every dimension
    g_total = high(box) - c - g_sum
    @inbounds for i in 1:n
        if !iszero(g_total[i])
            g = SingleEntryVector(i, n, g_total[i])
            push!(G, g)
        end
    end

    return Zonotope(c, G)
end

"""
    overapproximate(lm::LinearMap{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate a lazy linear map of a zonotopic set with a zonotope.

### Input

- `lm`       -- lazy linear map of a zonotopic set
- `Zonotope` -- target set type

### Output

The tight zonotope corresponding to `lm`.
"""
function overapproximate(lm::LinearMap{N, <:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    return convert(Zonotope, lm)
end

# Given center, (dependent) generator matrix and exponent matrix of a (simple)
# sparse polynomial zonotope, compute thew new center and generator matrix of
# its zonotope overapproximation.
function _zonotope_overapprox!(c, G, E)
    @inbounds for (j, g) in enumerate(eachcol(G))
        if all(iseven, E[:, j])
            c .+= g / 2
            G[:, j] ./= 2
        end
    end
    return c, G
end

function _zonotope_overapprox(c, G, E)
    return _zonotope_overapprox!(copy(c), copy(G), E)
end

"""
    overapproximate(P::SimpleSparsePolynomialZonotope, ::Type{<:Zonotope};
                    [nsdiv]=1, [partition]=nothing)

Overapproximate a simple sparse polynomial zonotope with a zonotope.

### Input

- `P`         -- simple sparse polynomial zonotope
- `Zonotope`  -- target set type
- `nsdiv`     -- (optional, default: `1`) size of uniform partitioning grid
- `partition` -- (optional, default: `nothing`) tuple of integers indicating the
                 number of blocks in each dimensino; the length should match
                 `nparams(P)`

### Output

A zonotope.
"""
function overapproximate(P::SimpleSparsePolynomialZonotope, ::Type{<:Zonotope};
                         nsdiv=1, partition=nothing)
    if !isnothing(partition) || nsdiv != 1
        return overapproximate(P, UnionSetArray{Zonotope}; nsdiv=nsdiv,
                               partition=partition)
    end

    cnew, Gnew = _zonotope_overapprox(center(P), genmat(P), expmat(P))
    return Zonotope(cnew, Gnew)
end

"""
    overapproximate(P::SimpleSparsePolynomialZonotope, ::Type{<:Zonotope},
                    dom::IntervalBox)

Overapproximate a simple sparse polynomial zonotope over the parameter domain
`dom` with a zonotope.

### Input

- `P`         -- simple sparse polynomial zonotope
- `Zonotope`  -- target set type
- `dom`       -- parameter domain, which should be a subset of `[-1, 1]^q`,
                 where `q = nparams(P)`

### Output

A zonotope.
"""
function overapproximate(P::SimpleSparsePolynomialZonotope, ::Type{<:Zonotope},
                         dom::IA.IntervalBox)
    @assert dom ⊆ IA.IntervalBox(IA.Interval(-1, 1), nparams(P)) "dom should " *
        "be a subset of [-1, 1]^q"

    G = genmat(P)
    E = expmat(P)
    cnew = copy(center(P))
    Gnew = similar(G)
    @inbounds for (j, g) in enumerate(eachcol(G))
        # monomial value over the domain
        # α = mapreduce(x -> _fast_interval_pow(x[1],  x[2]), *, zip(dom, E[:, i]))
        α = IA.Interval(1, 1)
        for (i, vi) in enumerate(dom)
            α *= fast_interval_pow(vi, E[i, j])
        end
        m, r = IA.midpoint_radius(α)
        cnew .+= m * g
        Gnew[:, j] .= abs.(g) * r
    end
    return Zonotope(cnew, Gnew)
end

"""
    overapproximate(P::SparsePolynomialZonotope, ::Type{<:Zonotope})

Overapproximate a sparse polynomial zonotope with a zonotope.

### Input

- `P`         -- sparse polynomial zonotope
- `Zonotope`  -- target set type

### Output

A zonotope.
"""
function overapproximate(P::SparsePolynomialZonotope, ::Type{<:Zonotope})
    cnew, Gnew = _zonotope_overapprox(center(P), genmat_dep(P), expmat(P))
    return Zonotope(cnew, hcat(Gnew, genmat_indep(P)))
end

# function to be loaded by Requires
function load_taylormodels_overapproximation()
return quote

using .TaylorModels: Taylor1, TaylorN, TaylorModel1, TaylorModelN,
                     polynomial, remainder, domain,
                     normalize_taylor, linear_polynomial,
                     constant_term, evaluate, mid, get_numvars,
                     HomogeneousPolynomial

@inline function get_linear_coeffs(p::Taylor1)
    if p.order == 0
        return zeros(eltype(p), 1)
    end
    return linear_polynomial(p).coeffs[2:2]
end

@inline function get_linear_coeffs(p::TaylorN)
    if p.order == 0
        n = get_numvars()
        return zeros(eltype(p), n)
    end
    return linear_polynomial(p).coeffs[2].coeffs
end

# compute the nonlinear part of the polynomial p, without truncation
@inline function _nonlinear_polynomial(p::Taylor1{T}) where {T}
    pnl = deepcopy(p)
    pnl.coeffs[1] = zero(T)
    if p.order > 0
        pnl.coeffs[2] = zero(T)
    end
    return pnl
end

@inline function _nonlinear_polynomial(p::TaylorN{T}) where {T}
    pnl = deepcopy(p)
    pnl.coeffs[1] = HomogeneousPolynomial([zero(T)])
    if p.order > 0
        pnl.coeffs[2] = HomogeneousPolynomial([zero(T)])
    end
    return pnl
end

"""
    overapproximate(vTM::Vector{TaylorModel1{T, S}}, ::Type{<:Zonotope};
                    [remove_zero_generators]::Bool=true
                    [normalize]::Bool=true) where {T, S}

Overapproximate a Taylor model in one variable with a zonotope.

### Input

- `vTM`       -- vector of `TaylorModel1`
- `Zonotope`  --  target set type
- `remove_zero_generators` -- (optional; default: `true`) flag to remove zero
                              generators of the resulting zonotope
- `normalize` -- (optional; default: `true`) flag to skip the normalization of
                 the Taylor models

### Output

A zonotope that overapproximates the range of the given Taylor model.

### Examples

If the polynomials are linear, this method exactly transforms to a zonotope.
The nonlinear case necessarily introduces overapproximation error.
Consider the linear case first:

```jldoctest oa_tm1
julia> using LazySets, TaylorModels

julia> const IA = IntervalArithmetic;

julia> I = IA.Interval(-0.5, 0.5) # interval remainder
[-0.5, 0.5]

julia> x₀ = IA.Interval(0.0) # expansion point
[0, 0]

julia> D = IA.Interval(-3.0, 1.0)
[-3, 1]

julia> p1 = Taylor1([2.0, 1.0], 2) # define a linear polynomial
 2.0 + 1.0 t + 𝒪(t³)

julia> p2 = Taylor1([0.9, 3.0], 2) # define another linear polynomial
 0.9 + 3.0 t + 𝒪(t³)

julia> vTM = [TaylorModel1(pi, I, x₀, D) for pi in [p1, p2]]
2-element Vector{TaylorModel1{Float64, Float64}}:
  2.0 + 1.0 t + [-0.5, 0.5]
  0.9 + 3.0 t + [-0.5, 0.5]
```

Here, `vTM` is a Taylor model vector, since each component is a Taylor model in
one variable (`TaylorModel1`). Using `overapproximate(vTM, Zonotope)` we can
compute its associated zonotope in generator representation:

```jldoctest oa_tm1
julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
2-element Vector{Float64}:
  1.0
 -2.1

julia> Matrix(genmat(Z))
2×3 Matrix{Float64}:
 2.0  0.5  0.0
 6.0  0.0  0.5
```

Note how the generators of this zonotope mainly consist of two pieces: one comes
from the linear part of the polynomials, and another one corresponds to the
interval remainder. This conversion gives the same upper and lower bounds as the
range evaluation using interval arithmetic:

```jldoctest oa_tm1
julia> X = box_approximation(Z)
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([1.0, -2.1], [2.5, 6.5])

julia> Y = evaluate(vTM[1], vTM[1].dom) × evaluate(vTM[2], vTM[2].dom)
[-1.5, 3.5] × [-8.60001, 4.40001]

julia> H = convert(Hyperrectangle, Y) # this IntervalBox is the same as X
Hyperrectangle{Float64, StaticArraysCore.SVector{2, Float64}, StaticArraysCore.SVector{2, Float64}}([1.0, -2.1000000000000005], [2.5, 6.500000000000001])
```
However, the zonotope returns better results if we want to approximate the
Taylor model because it is not axis-aligned:

```jldoctest oa_tm1
julia> d = [-0.35, 0.93];

julia> ρ(d, Z) < ρ(d, X)
true
```

This method also works if the polynomials are non-linear; for example suppose
that we add a third polynomial with a quadratic term:

```jldoctest oa_tm1
julia> p3 = Taylor1([0.9, 3.0, 1.0], 3)
 0.9 + 3.0 t + 1.0 t² + 𝒪(t⁴)

julia> vTM = [TaylorModel1(pi, I, x₀, D) for pi in [p1, p2, p3]]
3-element Vector{TaylorModel1{Float64, Float64}}:
           2.0 + 1.0 t + [-0.5, 0.5]
           0.9 + 3.0 t + [-0.5, 0.5]
  0.9 + 3.0 t + 1.0 t² + [-0.5, 0.5]

julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
3-element Vector{Float64}:
  1.0
 -2.1
  2.4

julia> Matrix(genmat(Z))
3×4 Matrix{Float64}:
 2.0  0.5  0.0  0.0
 6.0  0.0  0.5  0.0
 6.0  0.0  0.0  5.0
```

The last generator corresponds to the addition of the interval remainder and the
box overapproximation of the nonlinear part of `p3` over the domain.

### Algorithm

Let ``\\text{vTM} = (p, I)`` be a vector of ``m`` Taylor models, where ``I``
is the interval remainder in ``\\mathbb{R}^m``. Let ``p_{lin}``
(resp. ``p_{nonlin}``) correspond to the linear (resp. nonlinear) part of each
scalar polynomial.

The range of ``\\text{vTM}`` can be enclosed by a zonotope with center ``c``
and matrix of generators ``G``, ``Z = ⟨c, G⟩``, by performing a conservative
linearization of ``\\text{vTM}``:

```math
    vTM' = (p', I') := (p_{lin} − p_{nonlin} , I + \\text{Int}(p_{nonlin})).
```

This algorithm proceeds in two steps:

1- Conservatively linearize ``\\text{vTM}`` as above and compute a box
   overapproximation of the nonlinear part.
2- Transform the linear Taylor model to a zonotope exactly through variable
   normalization onto the symmetric intervals ``[-1, 1]``.
"""
function overapproximate(vTM::Vector{TaylorModel1{T, S}}, ::Type{<:Zonotope};
                         remove_zero_generators::Bool=true,
                         normalize::Bool=true) where {T, S}
    return _overapproximate_vTM_zonotope(vTM, 1, T;
                                         remove_zero_generators=remove_zero_generators,
                                         normalize=normalize)
end

"""
    overapproximate(vTM::Vector{TaylorModelN{N, T, S}}, ::Type{<:Zonotope};
                    [remove_zero_generators]::Bool=true
                    [normalize]::Bool=true) where {N, T, S}


Overapproximate a multivariate Taylor model with a zonotope.

### Input

- `vTM`       -- vector of `TaylorModelN`
- `Zonotope`  -- target set type
- `remove_zero_generators` -- (optional; default: `true`) flag to remove zero
                              generators of the resulting zonotope
- `normalize` -- (optional; default: `true`) flag to skip the normalization of
                 the Taylor models

### Output

A zonotope that overapproximates the range of the given Taylor model.

### Examples

Consider a vector of two 2-dimensional Taylor models of order 2 and 4,
respectively.

```jldoctest
julia> using LazySets, TaylorModels

julia> const IA = IntervalArithmetic;

julia> x₁, x₂ = set_variables(Float64, ["x₁", "x₂"], order=8)
2-element Vector{TaylorN{Float64}}:
  1.0 x₁ + 𝒪(‖x‖⁹)
  1.0 x₂ + 𝒪(‖x‖⁹)

julia> x₀ = IntervalBox(0..0, 2) # expansion point
[0, 0]²

julia> Dx₁ = IA.Interval(0.0, 3.0) # domain for x₁
[0, 3]

julia> Dx₂ = IA.Interval(-1.0, 1.0) # domain for x₂
[-1, 1]

julia> D = Dx₁ × Dx₂ # take the Cartesian product of the domain on each variable
[0, 3] × [-1, 1]

julia> r = IA.Interval(-0.5, 0.5) # interval remainder
[-0.5, 0.5]

julia> p1 = 1 + x₁^2 - x₂
 1.0 - 1.0 x₂ + 1.0 x₁² + 𝒪(‖x‖⁹)

julia> p2 = x₂^3 + 3x₁^4 + x₁ + 1
 1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + 𝒪(‖x‖⁹)

julia> vTM = [TaylorModelN(pi, r, x₀, D) for pi in [p1, p2]]
2-element Vector{TaylorModelN{2, Float64, Float64}}:
            1.0 - 1.0 x₂ + 1.0 x₁² + [-0.5, 0.5]
  1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + [-0.5, 0.5]

julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
2-element Vector{Float64}:
   5.5
 124.0

julia> Matrix(genmat(Z))
2×4 Matrix{Float64}:
 0.0  -1.0  5.0    0.0
 1.5   0.0  0.0  123.0
```

### Algorithm

We refer to the algorithm description for the univariate case.
"""
function overapproximate(vTM::Vector{TaylorModelN{N, T, S}}, ::Type{<:Zonotope};
                         remove_zero_generators::Bool=true,
                         normalize::Bool=true) where {N, T, S}
    n = N  # number of variables is get_numvars() in TaylorSeries
    return _overapproximate_vTM_zonotope(vTM, n, T;
                                         remove_zero_generators=remove_zero_generators,
                                         normalize=normalize)
end

function _overapproximate_vTM_zonotope(vTM, n, N;
                                       remove_zero_generators::Bool=true,
                                       normalize::Bool=true)
    m = length(vTM)

    # preallocations
    c = Vector{N}(undef, m)  # center of the zonotope
    G = Matrix{N}(undef, m, n + m)  # generator matrix

    @inbounds for (i, p) in enumerate(vTM)
        pol, dom = polynomial(p), domain(p)

        # linearize the TM
        pol_lin = constant_term(pol) + linear_polynomial(pol)

        # normalize the linear polynomial to the symmetric interval [-1, 1]^m
        pol_lin_norm = normalize ? normalize_taylor(pol_lin, dom, true) : pol_lin

        # overapproximate the nonlinear terms with an interval
        pol_nonlin = _nonlinear_polynomial(pol)
        rem_nonlin = evaluate(pol_nonlin, dom) + remainder(p)

        # build the generators
        α = mid(rem_nonlin)
        c[i] = constant_term(pol_lin_norm) + α  # constant terms
        G[i, 1:n] = get_linear_coeffs(pol_lin_norm)  # linear terms
        # interval generator
        for j in n+1:n+m
            G[i, j] = zero(N)
        end
        G[i, n+i] = abs(rem_nonlin.hi - α)
    end

    Z = Zonotope(c, G)
    if remove_zero_generators
        Z = LazySets.remove_zero_generators(Z)
    end
    return Z
end

end end  # quote / load_taylormodels_overapproximation

function load_intervalmatrices_overapproximation()
return quote

using .IntervalMatrices: AbstractIntervalMatrix, midpoint_radius

# temporary patch for IntervalArithmetic#317
function convert(::Type{IntervalMatrices.Interval{T}},
                 x::IntervalMatrices.Interval{T}) where {T<:Real}
    return x
end

"""
    overapproximate(lm::LinearMap{N, <:AbstractZonotope, NM,
                                  <:AbstractIntervalMatrix{NM}},
                    ::Type{<:Zonotope}) where {N, NM}

Overapproximate an interval-matrix linear map of a zonotopic set by a zonotope.

### Input

- `lm`       -- interval-matrix linear map of a zonotopic set
- `Zonotope` -- target set type

### Output

A zonotope overapproximating the linear map.

### Algorithm

This implementation uses the method proposed in [1].

Given an interval matrix ``M = \\tilde{M} + ⟨-\\hat{M},\\hat{M}⟩`` (split into a
conventional matrix and a symmetric interval matrix) and a zonotope
``⟨c, g_1, …, g_m⟩``, we compute the resulting zonotope
``⟨\\tilde{M}c, \\tilde{M}g_1, …, \\tilde{M}g_m, v_1, …, v_n⟩`` where the
``v_j``, ``j = 1, …, n``, are defined as

```math
    v_j = \\begin{cases} 0 & i ≠ j \\\\
          \\hat{M}_j (|c| + \\sum_{k=1}^m |g_k|) & i = j. \\end{cases}
```

[1] Althoff, Stursberg, Buss. *Reachability analysis of linear systems with
uncertain parameters and inputs*. CDC 2007.
"""
function overapproximate(lm::LinearMap{N, <:AbstractZonotope, NM,
                                       <:AbstractIntervalMatrix{NM}},
                         ::Type{<:Zonotope}) where {N, NM}
    Mc, Ms = midpoint_radius(lm.M)
    Z = lm.X
    c = Mc * center(Z)
    n = dim(lm)
    nG = ngens(Z)
    G = zeros(N, n, nG + n)
    vector_sum = abs.(center(Z))
    @inbounds for (j, g) in enumerate(generators(Z))
        G[:, j] = Mc * g
        vector_sum += abs.(g)
    end
    @inbounds for i in 1:n
        row = @view Ms[i, :]
        G[i, i + nG] = dot(row, vector_sum)
    end
    return Zonotope(c, G)
end

end end  # quote / load_intervalmatrices_overapproximation()

"""
    overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                    dir::Union{AbstractDirections, Type{<:AbstractDirections}};
                    [algorithm]="vrep", kwargs...)

Overapproximate a set with a zonotope.

### Input

- `X`         -- set
- `Zonotope`  -- target set type
- `dir`       -- directions used for the generators
- `algorithm` -- (optional, default: `"vrep"`) algorithm used to compute the
                 overapproximation
- `kwargs`    -- further algorithm-specific options

### Output

A zonotope that overapproximates `X` and uses at most the generator directions
provided in `dir` (redundant directions will be ignored).

### Notes

Two algorithms are available:

- `"vrep"` -- Overapproximate a polytopic set with a zonotope of minimal total
  generator sum using only generators in the given directions. Under this
  constraint, the zonotope has the minimal sum of generator vectors. See the
  docstring of [`_overapproximate_zonotope_vrep`](@ref) for further details.

- `"cpa"` -- Overapproximate a polytopic set with a zonotope using a Cartesian
  decomposition into two-dimensional blocks. See the docstring of
  [`_overapproximate_zonotope_cpa`](@ref) for further details.
"""
function overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                         dir::Union{AbstractDirections, Type{<:AbstractDirections}};
                         algorithm="vrep", kwargs...)
    if algorithm == "vrep"
        return  _overapproximate_zonotope_vrep(X, dir, kwargs...)
    elseif algorithm == "cpa"
        cpa = _overapproximate_zonotope_cpa(X, dir)
        return convert(Zonotope, cpa)
    else
        throw(ArgumentError("algorithm $algorithm is not known"))
    end
end

# disambiguation
function overapproximate(∅::EmptySet, ::Type{<:Zonotope},
                         dir::AbstractDirections;
                         algorithm="vrep", kwargs...)
    return ∅
end
function overapproximate(∅::EmptySet, ::Type{<:Zonotope},
                         dir::Type{<:AbstractDirections};
                         algorithm="vrep", kwargs...)
    return ∅
end

"""
    _overapproximate_zonotope_vrep(X::LazySet{N},
                                   dir::AbstractDirections;
                                   solver=default_lp_solver(N)) where {N}

Overapproximate a polytopic set with a zonotope of minimal total generator sum
using only generators in the given directions.

### Input

- `X`      -- polytopic set
- `dir`    -- directions used for the generators
- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

### Output

A zonotope that overapproximates `X` and uses at most the directions provided in
`dir` (redundant directions will be ignored).
Under this constraint, the zonotope has the minimal sum of generator vectors.

### Notes

The algorithm only requires one representative of each generator direction and
their additive inverse (e.g., only one of `[1, 0]` and `[-1, 0]`) and assumes
that the directions are normalized.
We preprocess the directions in that respect.

### Algorithm

We solve a linear program parametric in the vertices ``v_j`` of `X` and the
directions ``d_k`` in `dir` presented in Section 4.2 in [1], adapting the
notation to the one used in this library.

```math
    \\min \\sum_{k=1}^l α_k \\
    s.t. \\
    c + \\sum_{k=1}^l b_{kj} * d_k = v_j \\quad \\forall j \\
    -α_k ≤ b_{kj} ≤ α_k \\quad \\forall k, j \\
    α_k ≥ 0 \\quad \\forall k
```

The resulting zonotope has center `c` and generators `α_k · d_k`.

Note that the first type of side constraints is vector-based and that the
nonnegativity constraints (last type) are not stated explicitly in [1].

[1] *Zonotopes as bounding volumes*. L. J. Guibas, A. T. Nguyen, L. Zhang.
    SODA 2003.
"""
function _overapproximate_zonotope_vrep(X::LazySet{N},
                                        dir::AbstractDirections;
                                        solver=default_lp_solver(N)) where {N}
    @assert isconvextype(typeof(X)) "this algorithm requires a convex " *
                                    "polytope as input"

    # TODO "normalization" here involves two steps: removing opposite directions
    # and normalizing the direction vector
    # for the latter we can use the normalization information from dispatch on
    # the directions (some typical such directions are normalized by definition)
    function normalize_directions(dir::AD) where {AD}
        dirs = Vector{Vector{N}}()
        for d in dir
            d_normalized = normalize(d)
            if (d_normalized ∉ dirs) && (-d_normalized ∉ dirs)
                push!(dirs, d_normalized)
            end
        end
        return dirs
    end

    n = dim(X)
    dirs = normalize_directions(dir)
    l = length(dirs)
    V = vertices_list(X)
    m = length(V)
    nvariables = l + n + l * m  # l 'α_k' + n 'p[i]' + lm 'b_lj'
    nconstraints = n * m + 2 * l * m  # nm 'p_j' + 2lm 'b_kj'

    obj = vcat(ones(l), zeros(nvariables - l))
    A = spzeros(N, nconstraints, nvariables)
    sense = vcat(fill('=', n * m), fill('<', 2 * l * m))
    b = zeros(N, nconstraints)
    r = 1
    # constraints p[i] + \sum_k v_j*b_kj = p_j
    col_offset_p = l
    col_offset_b = l + n + 1
    for (j, vj) in enumerate(V)
        for i in 1:n
            A[r, col_offset_p + i] = N(1)  # p[i]
            for (k, dk) in enumerate(dirs)
                A[r, col_offset_b + (k - 1) * m] = dk[i]
            end
            b[r] = vj[i]
            r += 1
        end
        col_offset_b += 1
    end
    # constraints -α_k ± b_kj <= 0
    col_offset_b = l + n
    for k in 1:l
        for j in 1:m
            for sign_b in N[1, -1]
                A[r, k] = N(-1)  # - α_k
                A[r, col_offset_b + j] = sign_b  # ± b_kj
                r += 1
            end
        end
        col_offset_b += m
    end
    lbounds = vcat(zeros(N, l), fill(N(-Inf), nvariables - l))
    ubounds = fill(Inf, nvariables)
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)

    if !is_lp_optimal(lp.status)
        error("got unexpected status from LP solver: $(lp.status)")
    end
    c = lp.sol[(col_offset_p+1):(col_offset_p+n)]
    ck = lp.sol[1:l]
    G = Matrix{N}(undef, n, l)
    for (j, dj) in enumerate(dirs)
        G[:, j] = ck[j] * dj
    end
    return Zonotope(c, G)
end

# overload on direction type
function _overapproximate_zonotope_vrep(X::LazySet{N},
                                        dir::Type{<:AbstractDirections};
                                        solver=default_lp_solver(N)) where {N}
    return _overapproximate_zonotope_vrep(X, _get_directions(dir, dim(X)), solver=solver)
end

"""
    _overapproximate_zonotope_cpa(X::LazySet, dir::Type{<:AbstractDirections})

Overapproximate a polytopic set with a zonotope using Cartesian decomposition.

### Input

- `X`   -- polytopic set
- `dir` -- directions used for the generators

### Output

A zonotope that overapproximates `X`.

### Notes

The algorithm decomposes `X` into 2D sets and overapproximates those sets with
zonotopes, and finally converts the Cartesian product of the sets to a zonotope.

### Algorithm

The implementation is based on Section 8.2.4 in [1].

[1] Le Guernic, C. *Reachability analysis of hybrid systems with linear
continuous dynamics* (Doctoral dissertation). 2009.
"""
function _overapproximate_zonotope_cpa(X::LazySet, dir::AbstractDirections)
    n = dim(X)
    if dim(dir) != 2
        # try to convert to 2D directions
        dir = _get_directions(typeof(dir), 2)
    end
    # overapproximate 2D blocks
    if n > 1
        πX_2D = [project(X, [i, i+1]) for i in 1:2:n-1]
        Z_2D = [_overapproximate_zonotope_vrep(poly, dir) for poly in πX_2D]
    end

    if iseven(n)
        out = Z_2D
    else
        # odd case projects onto an interval
        πX_n = project(X, [n])
        πX_n = overapproximate(πX_n, Interval)
        πX_n = convert(Zonotope, πX_n)

        if n == 1
            out = [πX_n]
        else
            out = vcat(Z_2D, πX_n)
        end
    end
    return CartesianProductArray(out)
end

# overload on direction type
function _overapproximate_zonotope_cpa(X::LazySet, dir::Type{<:AbstractDirections})
    _overapproximate_zonotope_cpa(X, _get_directions(dir, 2))
end

"""
    overapproximate(r::Rectification{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the rectification of a zonotopic set with a zonotope.

### Input

- `r`        -- lazy rectification of a zonotopic set
- `Zonotope` -- target set type

### Output

A zonotope overapproximation of the set obtained by rectifying `Z`.

### Algorithm

This method implements [Theorem 3.1, 1].

[1] Singh, G., Gehr, T., Mirman, M., Püschel, M., & Vechev, M. *Fast and
effective robustness certification*. NeurIPS 2018.
"""
function overapproximate(r::Rectification{N, <:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    Z = set(r)
    c = copy(center(Z))
    G = copy(genmat(Z))
    n, m = size(G)
    row_idx = Vector{Int}()
    μ_idx = Vector{N}()

    @inbounds for i in 1:n
        lx, ux = low(Z, i), high(Z, i)
        if !_leq(lx, zero(N))
            nothing
        elseif _leq(ux, zero(N))
            c[i] = zero(N)
            for j in 1:m
                G[i, j] = zero(N)
            end
        else
            λ = ux / (ux - lx)
            μ = - λ * lx / 2
            c[i] = c[i] * λ + μ
            for j in 1:m
                G[i, j] = G[i, j] * λ
            end
            push!(row_idx, i)
            push!(μ_idx, μ)
        end
    end

    q = length(row_idx)
    if q >= 1
        Gnew = zeros(N, n, q)
        j = 1
        @inbounds for i in row_idx
            Gnew[i, j] = μ_idx[j]
            j += 1
        end
        Gout = hcat(G, Gnew)
    else
        Gout = G
    end

    return Zonotope(c, remove_zero_columns(Gout))
end

"""
    overapproximate(CHA::ConvexHullArray{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the convex hull of a list of zonotopic sets with a zonotope.

### Input

- `CHA`      -- convex hull array of zonotopic sets
- `Zonotope` -- target set type

### Output

A zonotope overapproximation of the convex hull array of zonotopic sets.

### Algorithm

This method iteratively applies the overapproximation algorithm to the
convex hull of two zonotopic sets from the given array of zonotopic sets.
"""
function overapproximate(CHA::ConvexHullArray{N, <:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    arr = array(CHA)
    n = length(arr)
    if n == 1
        return convert(Zonotope, arr[1])
    else
        Zaux = overapproximate(ConvexHull(arr[1], arr[2]), Zonotope)
        for k in 3:n
            Zaux = overapproximate(ConvexHull(Zaux, arr[k]), Zonotope)
        end
        return Zaux
    end
end

"""
    overapproximate(QM::QuadraticMap{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate a quadratic map of a zonotopic set with a zonotope.

### Input

- `QM`       -- quadratic map of a zonotopic set
- `Zonotope` -- target set type

### Output

A zonotope overapproximating the quadratic map of a zonotopic set.

### Notes

Mathematically, a quadratic map of a zonotope with matrices ``Q`` is defined as:

```math
Z_Q = \\right\\{ \\lambda | \\lambda_i = x^T Q\\^{(i)} x,~i = 1, \\ldots, n,~x \\in Z \\left\\}
```

### Algorithm

This method implements [Lemma 1, 1].

[1] Matthias Althoff and Bruce H. Krogh. *Avoiding geometric intersection
operations in reachability analysis of hybrid systems*. HSCC 2012.
"""
function overapproximate(QM::QuadraticMap{N, <:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    Z = QM.X
    G = genmat(Z)
    c = center(Z)
    n, p = size(G)
    h = Matrix{N}(undef, n, binomial(p+2, 2)-1)
    d = Vector{N}(undef, n)
    g(x) = view(G, :, x)
    cᵀ = c'
    for (i, Qᵢ) in enumerate(QM.Q)
        cᵀQᵢ = cᵀ * Qᵢ
        Qᵢc = Qᵢ * c
        aux = zero(N)
        for j=1:p
            aux += g(j)' * Qᵢ * g(j)
            h[i, j] = cᵀQᵢ * g(j) + g(j)' * Qᵢc
            h[i, p+j] = g(j)' * Qᵢ * g(j) / 2
        end
        d[i] = cᵀQᵢ * c + aux / 2
        l = 0
        for j=1:p-1
            gjᵀQᵢ = g(j)' * Qᵢ
            Qᵢgj = Qᵢ * g(j)
            for k=j+1:p
                l += 1
                h[i, 2p+l] = gjᵀQᵢ * g(k) + g(k)' * Qᵢgj
            end
        end
    end
    return Zonotope(d, remove_zero_columns(h))
end

"""
    overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                    ::Type{<:Zonotope})

Overapproximate the intersection of a zonotopic set and a hyperplane with a
zonotope.

### Input

- `X`        -- intersection of a zonotopic set and a hyperplane
- `Zonotope` -- target set type

### Output

A zonotope overapproximating the intersection.

### Algorithm

This method implements Algorithm 3 in [1].

[1] Moussa Maïga, Nacim Ramdani, Louise Travé-Massuyès, Christophe Combastel:
*A CSP versus a zonotope-based method for solving guard set intersection in
nonlinear hybrid reachability*. Mathematics in Computer Science (8) 2014.
"""
function overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                         ::Type{<:Zonotope}) where {N}
    return _overapproximate_zonotope_hyperplane(X.X, X.Y)
end

# symmetric method
function overapproximate(X::Intersection{N, <:Hyperplane, <:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    return _overapproximate_zonotope_hyperplane(X.Y, X.X)
end

function _overapproximate_zonotope_hyperplane(Z::AbstractZonotope, H::Hyperplane)
    c, G = center(Z), genmat(Z)
    a, b = H.a, H.b

    s = G' * a
    d = b - dot(a, c)
    @static if VERSION >= v"1.7"
        sT = s'
    else
        sT = s' .+ 0  # convert to Vector (`nullspace` fails for lazy transpose)
    end
    V0 = nullspace(sT)

    cs = s * d / dot(s, s)
    Gs = V0 * V0'

    c_cap = c + G * cs
    G_cap = G * Gs
    return Zonotope(c_cap, G_cap)
end
