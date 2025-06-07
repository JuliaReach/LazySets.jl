```@meta
CurrentModule = LazySets.ZonotopeModule
```

# [ZonotopeMD](@id def_ZonotopeMD)

```@docs
ZonotopeMD
```

## Conversion

```julia
Zonotope(Z::ZonotopeMD)
```

## Operations

```@docs
genmat(::ZonotopeMD)
cartesian_product(::ZonotopeMD, ::ZonotopeMD)
```

Undocumented implementations:

* [`center`](@ref center)
* [`ngens`](@ref ngens)
