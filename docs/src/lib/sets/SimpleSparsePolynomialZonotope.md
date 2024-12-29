```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```

# [SimpleSparsePolynomialZonotope](@id def_SimpleSparsePolynomialZonotope)

```@docs
SimpleSparsePolynomialZonotope
```

```@meta
CurrentModule = LazySets
```

Alias:

```@docs
PolynomialZonotope
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
convex_hull(::LazySet)
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
convex_hull(::SimpleSparsePolynomialZonotope)
expmat(::SimpleSparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
ngens(::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
ngens(::SimpleSparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
rand(::Type{SimpleSparsePolynomialZonotope})
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
remove_redundant_generators(::AbstractZonotope)
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
remove_redundant_generators(::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{<:AbstractMatrix}, ::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{<:AbstractMatrix}, ::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
cartesian_product(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
cartesian_product(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets
```
```@docs; canonical=false
convex_hull(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
convex_hull(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_combination(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```
```@docs
linear_combination(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`center`](@ref center(::LazySet))
```@meta
CurrentModule = LazySets
```
* [`genmat`](@ref genmat(::AbstractZonotope))
```@meta
CurrentModule = LazySets.API
```
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
```@meta
CurrentModule = LazySets
```
* [`ngens_indep`](@ref ngens_indep(::AbstractPolynomialZonotope))
```@meta
CurrentModule = LazySets.API
```
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

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
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
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
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
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
* [`isconvextype`](@ref isconvextype(::Type{AbstractPolynomialZonotope}))
* [`order`](@ref order(::AbstractPolynomialZonotope))

Inherited from [`AbstractSparsePolynomialZonotope`](@ref):
* [`ngens_dep`](@ref ngens_dep(::AbstractSparsePolynomialZonotope))
* [`nparams`](@ref nparams(::AbstractSparsePolynomialZonotope))
