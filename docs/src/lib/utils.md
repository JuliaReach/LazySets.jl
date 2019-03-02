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
dot_zero
get_radius!
isinvertible
ispermutation
issquare
is_right_turn
is_tighter_same_dir_2D
nonzero_indices
samedir
sign_cadlag
_random_zero_sum_vector
remove_duplicates_sorted!
reseed
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
