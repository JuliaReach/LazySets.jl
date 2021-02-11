```@meta
CurrentModule = LazySets
```

# [Star](@id def_Star)

```@docs
Star
dim(::Star)
center(::Star)
predicate(::Star)
basis(::Star)
ρ(::AbstractVector, ::Star)
σ(::AbstractVector, ::Star)
an_element(X::Star)
isempty(X::Star)
isbounded(X::Star; cond_tol::Number=DEFAULT_COND_TOL)
∈(x::AbstractVector, S::Star)
vertices_list(X::Star; apply_convex_hull::Bool=true)
constraints_list(X::Star)
linear_map(M::AbstractMatrix, X::Star)
```
