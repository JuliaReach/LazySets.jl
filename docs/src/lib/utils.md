```@meta
CurrentModule = LazySets
```

# Utility functions

## Arrays module

```@docs
Arrays
cross_product(::AbstractMatrix{N}) where {N<:Real}
delete_zero_columns!
dot_zero
inner
is_cyclic_permutation
isinvertible
ispermutation
is_right_turn
issquare
nonzero_indices
remove_duplicates_sorted!
right_turn
samedir
SingleEntryVector
to_negative_vector
_up
_dr
_above
minmax
arg_minmax
```

## Functions and Macros

```@docs
an_element_helper
binary_search_constraints
get_radius!
is_tighter_same_dir_2D
sign_cadlag
_random_zero_sum_vector
rectify
require(::Symbol)
reseed
same_block_structure
substitute
substitute!
σ_helper
get_constrained_lowdimset
@neutral
@absorbing
@neutral_absorbing
@declare_array_version
@array_neutral
@array_absorbing
```

## Types

```@docs
CachedPair
StrictlyIncreasingIndices
```

## Inspection of set interfaces

```@docs
LazySets.subtypes(::Any, ::Bool)
implementing_sets
```

## Sampling

```@docs
_sample_unit_nsphere_muller!
_sample_unit_nball_muller!
LazySets._canonical_length(X::LazySet{N}) where {N<:Real}
sample
LazySets.Sampler{S<:LazySet, D<:Distribution}
LazySets._rejection_sampling!
```
