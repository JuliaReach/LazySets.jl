```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

# Utility functions

## Helpers for internal use only

### Functions and Macros

```@docs
an_element_helper
binary_search_constraints
cross_product(::AbstractMatrix{N}) where {N<:Real}
delete_zero_columns
dot_zero
get_radius!
inner
isinvertible
ispermutation
iscounterclockwise
issquare
is_right_turn
is_tighter_same_dir_2D
nonzero_indices
samedir
sign_cadlag
_random_zero_sum_vector
rectify
remove_duplicates_sorted!
require(::Symbol)
reseed
same_block_structure
substitute
substitute!
σ_helper
@neutral
@absorbing
@neutral_absorbing
@declare_array_version
@array_neutral
@array_absorbing
```

### Types

```@docs
CachedPair
Approximations.UnitVector
StrictlyIncreasingIndices
```

### Inspection of set interfaces

```@docs
LazySets.subtypes(::Any, ::Bool)
```
