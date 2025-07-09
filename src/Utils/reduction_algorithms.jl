"""
    AbstractReductionMethod

Abstract supertype for order-reduction methods of a zonotopic set.
"""
abstract type AbstractReductionMethod end

"""
    GIR05 <: AbstractReductionMethod

Zonotope order-reduction method from [Girard05](@citet).
"""
struct GIR05 <: AbstractReductionMethod end

"""
    COMB03 <: AbstractReductionMethod

Zonotope order-reduction method from [Combastel03](@citet).
"""
struct COMB03 <: AbstractReductionMethod end

"""
    ASB10 <: AbstractReductionMethod

Zonotope order-reduction method from [AlthoffSB10](@citet).
"""
struct ASB10 <: AbstractReductionMethod end

"""
    SRMB16 <: AbstractReductionMethod

Zonotope order-reduction method from [ScottRMB16](@citet).

### Fields

- `ϵ` -- (optional; default: `1e-6`) pivot threshold
- `δ` -- (optional; default: `1e-3`) volume threshold

### Notes

The method reorders the generator matrix using reduced row echelon form (rref) to the form
``[T ~ V]``, then iteratively removes one generator from ``V`` while updating ``T``.
"""
struct SRMB16{N<:Number} <: AbstractReductionMethod
    ϵ::N
    δ::N

    SRMB16(ϵ::N=1e-6, δ::N=1e-3) where {N<:Number} = new{N}(ϵ, δ)
end
