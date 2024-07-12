```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```

# [SimpleSparsePolynomialZonotope](@id def_SimpleSparsePolynomialZonotope)

```@docs
SimpleSparsePolynomialZonotope
```

## Operations

```@docs
center(::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope)
expmat(::SimpleSparsePolynomialZonotope)
genmat(::SimpleSparsePolynomialZonotope)
ngens(::SimpleSparsePolynomialZonotope)
rand(::Type{SimpleSparsePolynomialZonotope})
remove_redundant_generators(::SimpleSparsePolynomialZonotope)
linear_map(::AbstractMatrix, ::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
```

```@meta
CurrentModule = LazySets
```

Alias:

```@docs
PolynomialZonotope
```

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`order`](@ref order(::AbstractPolynomialZonotope))

Inherited from [`AbstractSparsePolynomialZonotope`](@ref):
* [`ngens_dep`](@ref ngens_dep(::AbstractPolynomialZonotope))
* [`ngens_indep`](@ref ngens_indep(::AbstractPolynomialZonotope))
* [`nparams`](@ref nparams(::AbstractPolynomialZonotope))
