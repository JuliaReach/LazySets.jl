```@meta
CurrentModule = LazySets
```

# Utility functions

## Arrays module

```@docs
Arrays
cross_product(::AbstractMatrix{N}) where {N<:Real}
nonzero_columns
dot_zero
hasfullrowrank
inner
is_cyclic_permutation
isinvertible
ispermutation
is_right_turn
issquare
nonzero_indices
remove_duplicates_sorted!
remove_zero_columns
right_turn
samedir
SingleEntryVector
to_negative_vector
_up
_dr
_above
minmax
arg_minmax
extend
projection_matrix
LazySets.Arrays._vector_type
LazySets.Arrays._matrix_type
LazySets.Arrays.distance(::AbstractVector, ::AbstractVector, ::Real=2.0)
```

## Functions and Macros

```@docs
an_element_helper
binary_search_constraints
get_radius!
is_tighter_same_dir_2D
sign_cadlag
_leq_trig
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
sample
LazySets.Sampler
LazySets.RejectionSampler{S<:LazySet, D<:Distribution}
LazySets._sample!
```

## Volume

```@docs
volume
```
