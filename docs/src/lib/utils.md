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
ismultiple
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
extend
projection_matrix
LazySets.Arrays._vector_type
LazySets.Arrays._matrix_type
LazySets.Arrays.allequal
LazySets.Arrays.distance(::AbstractVector, ::AbstractVector; ::Real=2.0)
LazySets.Arrays.ispermutation
LazySets.Arrays.same_sign
LazySets.Arrays.to_matrix
LazySets.Arrays._rationalize
LazySets.Arrays.substitute
LazySets.Arrays.substitute!
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
same_block_structure
σ_helper
get_constrained_lowdimset
minmax
arg_minmax
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
```

## Inspection of set interfaces

```@docs
LazySets.subtypes(::Any, ::Bool)
implementing_sets
```

## Reading and writing

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
```

## Volume

```@docs
volume
```

## Symbolics

```@docs
LazySets._vec
```
