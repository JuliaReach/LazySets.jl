```@meta
CurrentModule = LazySets.SimpleSparsePolynomialZonotopeModule
```

# [SimpleSparsePolynomialZonotope](@id def_SimpleSparsePolynomialZonotope)

```@docs
SimpleSparsePolynomialZonotope
```

## Conversion

```julia
convert(::Type{SimpleSparsePolynomialZonotope}, ::AbstractSparsePolynomialZonotope)
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
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))

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
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isconvextype`](@ref isconvextype(::Type{AbstractPolynomialZonotope}))
* [`isoperation`](@ref isoperation(::LazySet))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`ispolyhedraltype`](@ref ispolyhedraltype(::Type{LazySet}))
* [`ispolytopic`](@ref ispolytopic(::LazySet))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`issubset`](@ref issubset(::LazySet, ::LazySet))
* [`linear_combination`](@ref linear_combination(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`center`](@ref center(::AbstractPolynomialZonotope, ::Int))
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`extrema`](@ref extrema(::AbstractPolynomialZonotope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolynomialZonotope}))
* [`isempty`](@ref isempty(::AbstractPolynomialZonotope))
* [`isuniversal`](@ref isuniversal(::AbstractPolynomialZonotope))
* [`order`](@ref order(::AbstractPolynomialZonotope))

Inherited from [`AbstractSparsePolynomialZonotope`](@ref):
* [`ngens_dep`](@ref ngens_dep(::AbstractSparsePolynomialZonotope))
* [`nparams`](@ref nparams(::AbstractSparsePolynomialZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractSparsePolynomialZonotope))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractSparsePolynomialZonotope, ::AbstractSparsePolynomialZonotope))
