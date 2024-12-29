```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```

# [SparsePolynomialZonotope](@id def_SparsePolynomialZonotope)

```@docs
SparsePolynomialZonotope
```

## Operations

```@docs
expmat(::SparsePolynomialZonotope)
genmat_dep(::SparsePolynomialZonotope)
genmat_indep(::SparsePolynomialZonotope)
indexvector(::SparsePolynomialZonotope)
polynomial_order(::SparsePolynomialZonotope)
rand(::Type{SparsePolynomialZonotope})
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
remove_redundant_generators(::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```
```@docs
remove_redundant_generators(::SparsePolynomialZonotope)
uniqueID(::Int)
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
reduce_order(::AbstractZonotope, ::Real)
```
```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```
```@docs
reduce_order(::SparsePolynomialZonotope, ::Real, ::AbstractReductionMethod=GIR05())
ρ(::AbstractVector, ::SparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
cartesian_product(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```
```@docs
cartesian_product(::SparsePolynomialZonotope, ::SparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
exact_sum(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```
```@docs
exact_sum(::SparsePolynomialZonotope, ::SparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_sum(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.SparsePolynomialZonotopeModule
```
```@docs
minkowski_sum(::SparsePolynomialZonotope, ::SparsePolynomialZonotope)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`center`](@ref center(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`an_element`](@ref an_element(::LazySet))
* [`area`](@ref area(::LazySet))
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`σ`](@ref σ(::AbstractVector, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`linear_combination`](@ref linear_combination(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolynomialZonotope}))
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`ngens`](@ref ngens(::AbstractPolynomialZonotope))
* [`order`](@ref order(::AbstractPolynomialZonotope))

Inherited from [`AbstractSparsePolynomialZonotope`](@ref):
* [`ngens_dep`](@ref ngens_dep(::AbstractSparsePolynomialZonotope))
* [`ngens_indep`](@ref ngens_indep(::AbstractSparsePolynomialZonotope))
* [`nparams`](@ref nparams(::AbstractSparsePolynomialZonotope))
