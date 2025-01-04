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
```

Undocumented implementations:

* [`ngens_dep`](@ref ngens_dep(::AbstractPolynomialZonotope))
* [`ngens_indep`](@ref ngens_indep(::AbstractPolynomialZonotope))

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
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
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
