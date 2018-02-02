export AbstractNonConvexSet,
       PolynomialZonotope,
       dim,
       polynomial_order,
       order,
       linear_map,
       scale,
       minkowsi_sum

"""
    AbstractNonConvexSet

Abstract type for non-convex sets.
"""
abstract type AbstractNonConvexSet{N} end

"""
    PolynomialZonotope{N} <: AbstractNonConvexSet{N}

Type that represents a polynomial zonotope.

### Fields

- `c`  -- starting point
- `E`  -- multi-indexed generators such that *all* indices have the same value
- `F`  -- multi-indexed generators such that *not all* indices have the same value
- `G`  -- single-indexed generators

### Notes

Polynomial zonotopes were introduced by M. Althoff in *Reachability analysis of nonlinear
systems using conservative polynomialization and non-convex sets*, Hybrid Systems:
Computation and Control, 2013, pp. 173–182.

### Examples

```jldoctest
julia> c = zeros(2);
julia> E1, E2 = diagm([-1, 0.5]), [1 1; 0.5 0.3];
julia> E = [E1, E2];
julia> F2 = [-0.5 1]';
julia> F = [F2];
julia> G = diagm([0.3, 0.3]);
julia> p = PolynomialZonotope(c, E, F, G)
LazySets.PolynomialZonotope{Float64}([0.0, 0.0], Array{Float64,2}[[-1.0 0.0; 0.0 0.5], [1.0 1.0; 0.5 0.3]], Array{Float64,2}[[-0.5; 1.0]], [0.3 0.0; 0.0 0.3])
julia> dim(p)
2
julia> order(p)
7//2
julia> polynomial_order(p)
2
```
"""
struct PolynomialZonotope{N} <: AbstractNonConvexSet{N}
    c::Vector{N}
    E::Vector{Matrix{N}}
    F::Vector{Matrix{N}}
    G::Matrix{N}

    # default constructor with dimension check
    function PolynomialZonotope{N}(c::Vector{N},
                                   E::Vector{Matrix{N}},
                                   F::Vector{Matrix{N}},
                                   G::Matrix{N}) where {N}

        # check polynomial order
        @assert length(E) == 1 + length(F)

        return new(c, E, F, G)
    end
end

# type-less convenience constructor
PolynomialZonotope(c::Vector{N},
                   E::Vector{Matrix{N}},
                   F::Vector{Matrix{N}},
                   G::Matrix{N}) where {N} = PolynomialZonotope{N}(c, E, F, G)

"""
    dim(pz::PolynomialZonotope)::Int

Return the ambient dimension of a polynomial zonotope.

### Input

- `pz` -- polynomial zonotope

### Output

An integer representing the ambient dimension of the polynomial zonotope.
"""
dim(pz::PolynomialZonotope)::Int = length(pz.c)

"""
    polynomial_order(pz::PolynomialZonotope)::Int

Polynomial order of a polynomial zonotope.

### Input

- `pz` -- polynomial zonotope

## Output

The polynomial order, defined as the maximal power of the scalar factors ``β_i``.
Usually denoted ``η``.
"""
polynomial_order(pz::PolynomialZonotope)::Int = length(pz.E)

"""
    order(pz::PolynomialZonotope)::Rational{Int}

Order of a polynomial zonotope.

### Input

- `pz` -- polynomial zonotope

## Output

The order, defined as the number of generators divided by the ambient dimension.
"""
function order(pz::PolynomialZonotope)::Rational{Int}
    η = polynomial_order(pz)  # number of generators
    p = size(pz.E[1], 2)  # number of dependent factors
    q = size(pz.G, 2)  # number of independent factors
    ξ = sum(i -> binomial(p+i-1, i), 1:η) + q
    n = dim(pz)
    return ξ//n
end

"""
    linear_map(M::Matrix, pz::PolynomialZonotope)

Return the linear map of a polynomial zonotope.

### Input

- `M`  -- matrix
- `pz` -- polynomial zonotope

## Output

Polynomial zonotope such that its center and generators are those of `pz`
multiplied by the matrix `M`.
"""
function linear_map(M::Matrix, pz::PolynomialZonotope)
    c = M * pz.c
    E = [M*Ei for Ei in pz.E]
    F = [M*Ci for Fi in pz.F]
    G = M * pz.G
    return PolynomialZonotope(c, E, F, G)
end

"""
    scale(α::Number, pz::PolynomialZonotope)

Return a polynomial zonotope modified by a scalar factor.

### Input

- `α` -- polynomial zonotope
- `pz` -- polynomial zonotope

## Output

Polynomial zonotope such that its center and generators are multiples of those
of `pz` by a factor ``α``.
"""
function scale(α::Number, pz::PolynomialZonotope)
    c = α * pz.c
    E = [α*Ei for Ei in pz.E]
    F = [α*Ci for Fi in pz.F]
    G = α * pz.G
    return PolynomialZonotope(c, E, F, G)
end

"""
    minkowski_sum(pz::PolynomialZonotope, z::Zonotope)

Return the Minkowski sum of a polynomial zonotope and a usual zonotope.

### Input

- `pz` -- polynomial zonotope
- `z`  -- usual zonotope

## Output

Polynomial zonotope.
"""
function minkowski_sum(pz::PolynomialZonotope, z::Zonotope)
    c = pz.c + z.center
    G = [G z.generators]
    return PolynomialZonotope(c, pz.E, pz.F, G)
end
