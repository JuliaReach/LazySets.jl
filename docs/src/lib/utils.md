```@meta
CurrentModule = LazySets
```

# Utilities

## Macros

```@docs
@neutral
@absorbing
@neutral_absorbing
@declare_binary_operation
@declare_array_version
@array_neutral
@array_absorbing
```

## Types

```@docs
CachedPair
```

## Inspection of set interfaces

```@docs
implementing_sets
```

## File formats

```@docs
LazySets.read_gen(::String)
```

## Sampling

```@docs
_sample_unit_nsphere_muller!
_sample_unit_nball_muller!
sample
LazySets.AbstractSampler
LazySets.CombinedSampler
LazySets.FaceSampler
LazySets.HalfSpaceSampler
LazySets.HyperplaneSampler
LazySets.SingletonSampler
LazySets.RejectionSampler
LazySets.RandomWalkSampler
LazySets.PolynomialZonotopeSampler
```

## Symbolics

```@docs
LazySets._vec
```

## SymEngine

```@docs
LazySets.free_symbols
LazySets._is_linearcombination
LazySets._is_halfspace
LazySets._is_hyperplane
LazySets.convert(::Type{HalfSpace{N}}, ::Expr; vars::Vector{Basic}=Basic[]) where {N}
LazySets.convert(::Type{Hyperplane{N}}, ::Expr; vars::Vector{Basic}=Basic[]) where {N}
```

## Functions for numbers

```@docs
sign_cadlag
minmax
arg_minmax
```

## Other functions

```@docs
_an_element_helper_hyperplane
binary_search_constraints
get_constrained_lowdimset
get_radius!
is_tighter_same_dir_2D
_leq_trig
same_block_structure
_σ_hyperplane_halfspace
```
