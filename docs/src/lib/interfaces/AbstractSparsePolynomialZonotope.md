```@contents
Pages = ["AbstractSparsePolynomialZonotope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Sparse polynomial zonotope sets (AbstractSparsePolynomialZonotope)](@id def_AbstractSparsePolynomialZonotope)

```@docs
AbstractSparsePolynomialZonotope
```

This interface requires to implement the following functions:

```@docs
expmat(::AbstractSparsePolynomialZonotope)
genmat_dep(::AbstractSparsePolynomialZonotope)
genmat_indep(::AbstractSparsePolynomialZonotope)
```

This interface defines the following functions:

```@docs
nparams(::AbstractSparsePolynomialZonotope)
ρ(::AbstractVector, ::AbstractSparsePolynomialZonotope)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_combination(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
linear_combination(::AbstractSparsePolynomialZonotope, ::AbstractSparsePolynomialZonotope)
```

Undocumented implementations:

```@meta
CurrentModule = LazySets.API
```
* [`extrema`](@ref extrema(::LazySet))
```@meta
CurrentModule = LazySets
```
* [`ngens_dep`](@ref ngens_dep(::AbstractPolynomialZonotope))
* [`ngens_indep`](@ref ngens_indep(::AbstractPolynomialZonotope))
```@meta
CurrentModule = LazySets.API
```
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))

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
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`isbounded`](@ref isbounded(::LazySet))
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
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`triangulate`](@ref triangulate(::LazySet))
* [`triangulate_faces`](@ref triangulate_faces(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`isapprox`](@ref isapprox(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, ::LazySet))

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`center`](@ref center(::AbstractPolynomialZonotope, ::Int))
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolynomialZonotope}))
* [`isempty`](@ref isempty(::AbstractPolynomialZonotope))
* [`isuniversal`](@ref isuniversal(::AbstractPolynomialZonotope))
* [`ngens`](@ref ngens(::AbstractPolynomialZonotope))
* [`order`](@ref order(::AbstractPolynomialZonotope))

## Implementations

* [Sparse polynomial zonotope (SparsePolynomialZonotope)](@ref def_SparsePolynomialZonotope)
* [Simplified sparse polynomial zonotope (SimpleSparsePolynomialZonotope)](@ref def_SimpleSparsePolynomialZonotope)
