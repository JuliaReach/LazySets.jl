```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

# Utility functions

```@docs
sign_cadlag
ispermutation
@neutral
@absorbing
@declare_array_version
```

## Helpers for internal use only

### Functions and Macros

```@docs
@neutral_absorbing
@array_neutral
@array_absorbing
get_radius!
an_element_helper
σ_helper
binary_search_constraints
nonzero_indices
samedir
```

### Types

```@docs
CachedPair
Approximations.UnitVector
```
