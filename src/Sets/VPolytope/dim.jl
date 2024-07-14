"""
    dim(P::VPolytope)

Return the dimension of a polytope in vertex representation.

### Input

- `P` -- polytope in vertex representation

### Output

The ambient dimension of the polytope in vertex representation.
If `P` is empty, the result is ``-1``.

### Examples

```jldoctest
julia> P = VPolytope();

julia> isempty(P.vertices)
true

julia> dim(P)
-1

julia> P = VPolytope([ones(3)]);

julia> P.vertices
1-element Vector{Vector{Float64}}:
 [1.0, 1.0, 1.0]

julia> dim(P) == 3
true
```
"""
function dim(P::VPolytope)
    return isempty(P.vertices) ? -1 : @inbounds length(P.vertices[1])
end
