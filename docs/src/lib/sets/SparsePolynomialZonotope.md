```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```

# [SparsePolynomialZonotope](@id def_SparsePolynomialZonotope)

```@docs
SparsePolynomialZonotope
```

## Operations

```@docs
center(::SparsePolynomialZonotope)
expmat(::SparsePolynomialZonotope)
genmat_dep(::SparsePolynomialZonotope)
genmat_indep(::SparsePolynomialZonotope)
indexvector(::SparsePolynomialZonotope)
polynomial_order(::SparsePolynomialZonotope)
rand(::Type{SparsePolynomialZonotope})
remove_redundant_generators(::SparsePolynomialZonotope)
uniqueID(::Int)
linear_map(::AbstractMatrix, ::SparsePolynomialZonotope)
reduce_order(::SparsePolynomialZonotope, ::Real, ::AbstractReductionMethod=GIR05())
ρ(::AbstractVector, ::SparsePolynomialZonotope)
translate(::SparsePolynomialZonotope, ::AbstractVector)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`ngens`](@ref ngens(::AbstractPolynomialZonotope))
* [`order`](@ref order(::AbstractPolynomialZonotope))

Inherited from [`AbstractSparsePolynomialZonotope`](@ref):
* [`ngens_dep`](@ref ngens_dep(::AbstractPolynomialZonotope))
* [`ngens_indep`](@ref ngens_indep(::AbstractPolynomialZonotope))
* [`nparams`](@ref nparams(::AbstractPolynomialZonotope))
