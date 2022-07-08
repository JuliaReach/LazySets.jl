```@meta
CurrentModule = LazySets
```

# [Line](@id def_Line)

```@docs
Line
dim(::Line)
ρ(::AbstractVector, ::Line)
σ(::AbstractVector, ::Line)
∈(::AbstractVector, ::Line)
an_element(::Line)
direction(::Line)
rand(::Type{Line})
isbounded(::Line)
isuniversal(::Line; ::Bool=false)
isempty(::Line)
constraints_list(::Line{N, VN}) where {N, VN}
translate(::Line, ::AbstractVector)
translate!(::Line, ::AbstractVector)
normalize(::Line, ::Real=2.0)
normalize!(::Line, ::Real=2.0)
distance(::AbstractVector, ::Line; ::Real=2.0)
linear_map(::AbstractMatrix, ::Line)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
