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
get_radius!
ispermutation
nonzero_indices
samedir
sign_cadlag
_random_zero_sum_vector
remove_duplicates_sorted!
reseed
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
