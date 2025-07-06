```@meta
CurrentModule = LazySets
```

# Utilities

## Array set types

```@docs
flatten
neutral
absorbing
```

### Internal helper macros

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

## Polyhedra

```@docs
isfeasible(::AbstractMatrix, ::AbstractVector, ::Bool=false)
```

## File formats

```@docs
LazySets.read_gen(::String)
```

## Sampling

```@docs
_sample_unit_nsphere_muller!
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
LazySets.UniverseSampler
```

## Symbolics

```@docs
LazySets._vec
```

## SymEngine

```@docs
LazySets.free_symbols
LazySets._is_linear_combination
```

## Functions for numbers

```@docs
sign_cadlag
```

## Other functions

```@docs
binary_search_constraints
get_constrained_lowdimset
get_radius!
is_tighter_same_dir_2D
_leq_trig
same_block_structure
_σ_hyperplane_halfspace
```
