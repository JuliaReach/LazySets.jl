module LazySetsTaylorModelsExt

using IntervalArithmetic: Interval, inf, sup
import LazySets
using LazySets.HyperrectangleModule: Hyperrectangle
using LazySets.ZonotopeModule: Zonotope
using TaylorModels: HomogeneousPolynomial, TaylorModel1, TaylorModelN, Taylor1,
                    TaylorN, constant_term, domain, evaluate, get_numvars,
                    get_order, linear_polynomial, mid, normalize_taylor,
                    polynomial, remainder
using TaylorModels.TaylorSeries: numtype  # NOTE: this is an internal function
import LazySets.Approximations: box_approximation, overapproximate

include("TaylorModels/SparsePolynomialZonotopeExt.jl")

function _eltype_TM(TMi)
    return numtype(polynomial(TMi))
end

function box_approximation(vTM::Vector{<:TaylorModel1})
    return _box_approximation_vTM(vTM)
end

function box_approximation(vTM::Vector{<:TaylorModelN})
    return _box_approximation_vTM(vTM)
end

function _box_approximation_vTM(vTM)
    bounds = [evaluate(vTM[i], domain(p)) for (i, p) in enumerate(vTM)]
    return Hyperrectangle(; low=inf.(bounds), high=sup.(bounds))
end

# internal helper function
@inline function get_linear_coeffs(p::Taylor1)
    if get_order(p) == 0
        return zeros(eltype(p), 1)
    end
    return linear_polynomial(p).coeffs[2:2]
end

# internal helper function
@inline function get_linear_coeffs(p::TaylorN)
    if get_order(p) == 0
        n = get_numvars()
        return zeros(eltype(p), n)
    end
    return linear_polynomial(p).coeffs[2].coeffs
end

# internal helper function
# compute the nonlinear part of the polynomial p, without truncation
@inline function _nonlinear_polynomial(p::Taylor1{T}) where {T}
    pnl = deepcopy(p)
    pnl.coeffs[1] = zero(T)
    if get_order(p) > 0
        pnl.coeffs[2] = zero(T)
    end
    return pnl
end

# internal helper function
@inline function _nonlinear_polynomial(p::TaylorN{T}) where {T}
    pnl = deepcopy(p)
    pnl.coeffs[1] = HomogeneousPolynomial([zero(T)])
    if get_order(p) > 0
        pnl.coeffs[2] = HomogeneousPolynomial([zero(T)])
    end
    return pnl
end

"""
    overapproximate(vTM::Vector{TaylorModel1{T, S}}, ::Type{<:Zonotope};
                    [remove_redundant_generators]::Bool=true
                    [normalize]::Bool=true) where {T, S}

Overapproximate a Taylor model in one variable with a zonotope.

### Input

- `vTM`       -- vector of `TaylorModel1`
- `Zonotope`  --  target set type
- `remove_redundant_generators` -- (optional; default: `true`) flag to remove
                                    redundant generators of the resulting zonotope
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

julia> I = IA.interval(-0.5, 0.5) # interval remainder
[-0.5, 0.5]

julia> x₀ = IA.interval(0.0) # expansion point
[0.0, 0.0]

julia> D = IA.interval(-3.0, 1.0)
[-3.0, 1.0]

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

```@meta
DocTestSetup = quote
    # TODO needed while TaylorModels throws warnings; see also the corresponding block below
    import Logging
    Logging.disable_logging(Logging.Warn)
end
DocTestTeardown = quote
    Logging.disable_logging(Logging.Debug)
end
```

```jldoctest oa_tm1
julia> Z = overapproximate(vTM, Zonotope);

julia> center(Z)
2-element Vector{Float64}:
  1.0
 -2.0999999999999996

julia> Matrix(genmat(Z))
2×3 Matrix{Float64}:
 0.0  2.0  0.5
 0.5  6.0  0.0
```

Note how the generators of this zonotope mainly consist of two pieces: one comes
from the linear part of the polynomials, and another one corresponds to the
interval remainder. This conversion gives the same upper and lower bounds as the
range evaluation using interval arithmetic:

```jldoctest oa_tm1
julia> X = box_approximation(Z)
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([1.0, -2.0999999999999996], [2.5, 6.5])

julia> Y = [evaluate(vTM[1], vTM[1].dom), evaluate(vTM[2], vTM[2].dom)]
2-element Vector{IntervalArithmetic.Interval{Float64}}:
 [-1.5, 3.5]
 [-8.60001, 4.40001]
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
 -2.0999999999999996
  2.4000000000000004

julia> Matrix(genmat(Z))
3×4 Matrix{Float64}:
 0.0  0.0  2.0  0.5
 0.0  0.5  6.0  0.0
 5.0  0.0  6.0  0.0
```

```@meta
DocTestSetup = nothing
DocTestTeardown = nothing
```

The last generator corresponds to the addition of the interval remainder and the
box overapproximation of the nonlinear part of `p3` over the domain.

### Algorithm

Let ``\\text{vTM} = (p, I)`` be a vector of ``m`` Taylor models, where ``I``
is the interval remainder in ``ℝ^m``. Let ``p_{lin}``
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
function overapproximate(vTM::Vector{TaylorModel1{T,S}}, ::Type{<:Zonotope};
                            remove_redundant_generators::Bool=true,
                            normalize::Bool=true) where {T,S}
    return _overapproximate_vTM_zonotope(vTM, 1, T;
                                            remove_redundant_generators=remove_redundant_generators,
                                            normalize=normalize)
end

"""
    overapproximate(vTM::Vector{<:TaylorModelN{N, T}}, ::Type{<:Zonotope};
                    [remove_redundant_generators]::Bool=true
                    [normalize]::Bool=true) where {N,T}


Overapproximate a multivariate Taylor model with a zonotope.

### Input

- `vTM`       -- vector of `TaylorModelN`
- `Zonotope`  -- target set type
- `remove_redundant_generators` -- (optional; default: `true`) flag to remove
                                    redundant generators of the resulting zonotope
- `normalize` -- (optional; default: `true`) flag to skip the normalization of
                    the Taylor models

### Output

A zonotope that overapproximates the range of the given Taylor model.

### Examples

Consider a vector of two 2-dimensional Taylor models of order 2 and 4,
respectively.

```@meta
DocTestSetup = quote
    # TODO needed while TaylorModels throws warnings; see also the corresponding block below
    import Logging
    Logging.disable_logging(Logging.Warn)
end
DocTestTeardown = quote
    Logging.disable_logging(Logging.Debug)
end
```

```jldoctest
julia> using LazySets, TaylorModels

julia> const IA = IntervalArithmetic;

julia> x₁, x₂ = set_variables(Float64, ["x₁", "x₂"], order=8)
2-element Vector{TaylorN{Float64}}:
  1.0 x₁ + 𝒪(‖x‖⁹)
  1.0 x₂ + 𝒪(‖x‖⁹)

julia> x₀ = fill(0..0, 2) # expansion point
2-element Vector{IntervalArithmetic.Interval{Float64}}:
 [0.0, 0.0]
 [0.0, 0.0]

julia> Dx₁ = IA.interval(0.0, 3.0) # domain for x₁
[0.0, 3.0]

julia> Dx₂ = IA.interval(-1.0, 1.0) # domain for x₂
[-1.0, 1.0]

julia> D = [Dx₁, Dx₂]
2-element Vector{IntervalArithmetic.Interval{Float64}}:
  [0.0, 3.0]
 [-1.0, 1.0]

julia> r = IA.interval(-0.5, 0.5) # interval remainder
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
2×2 Matrix{Float64}:
   0.0   6.0
 124.5  -0.0
```

```@meta
DocTestSetup = nothing
DocTestTeardown = nothing
```

### Algorithm

We refer to the algorithm description for the univariate case.
"""
function overapproximate(vTM::Vector{<:TaylorModelN{N,T}}, ::Type{<:Zonotope};
                         remove_redundant_generators::Bool=true,
                         normalize::Bool=true) where {N,T}
    n = N  # number of variables is get_numvars() in TaylorSeries
    return _overapproximate_vTM_zonotope(vTM, n, T;
                                         remove_redundant_generators=remove_redundant_generators,
                                         normalize=normalize)
end

function _isthin_approx(x::Interval)
    return inf(x) ≈ sup(x)
end

function _overapproximate_vTM_zonotope(vTM, n, N;
                                       remove_redundant_generators::Bool=true,
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
        cterm = constant_term(pol_lin_norm)
        @assert _isthin_approx(cterm) "unexpected interval value: $cterm"
        c[i] = mid(cterm) + α  # constant terms
        lin_coeffs = get_linear_coeffs(pol_lin_norm)
        @assert all(_isthin_approx, lin_coeffs) "unexpected interval value: $lin_coeffs"
        G[i, 1:n] = mid.(lin_coeffs)  # linear terms
        # interval generator
        for j in (n + 1):(n + m)
            G[i, j] = zero(N)
        end
        G[i, n + i] = abs(sup(rem_nonlin) - α)
    end

    if remove_redundant_generators
        G = LazySets.remove_redundant_generators(G)
    end
    return Zonotope(c, G)
end

end  # module
