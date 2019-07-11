```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets, Distributions
end
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
issquare
nonzero_indices
remove_duplicates_sorted!
samedir
SingleEntryVector
to_negative_vector
_up
_dr
_above
```

## Functions and Macros

```@docs
an_element_helper
binary_search_constraints
get_radius!
is_right_turn
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
```

## Sampling

```@docs
LazySets._sample_unit_nsphere_muller!
LazySets._sample_unit_nball_muller!
```
