"""
# Extended help

    ∈(x::AbstractVector, P::VPolygon)

### Algorithm

This implementation exploits that the polygon's vertices are sorted in
counter-clockwise fashion.
Under this assumption we can just check if the vertex lies on the left of each
edge, using the dot product.

### Examples

```jldoctest
julia> P = VPolygon([[2.0, 3.0], [3.0, 1.0], [5.0, 1.0], [4.0, 5.0]]);

julia> [4.5, 3.1] ∈ P
false
julia> [4.5, 3.0] ∈ P
true
julia> [4.4, 3.4] ∈ P  #  point lies on the edge
true
```
"""
@validate function ∈(x::AbstractVector, P::VPolygon)
    # special cases: 0 or 1 vertex
    @inbounds begin
        if length(P.vertices) == 0
            return false
        elseif length(P.vertices) == 1
            return x == P.vertices[1]
        end

        if !is_right_turn(P.vertices[1], x, P.vertices[end])
            return false
        end
        for i in 2:length(P.vertices)
            if !is_right_turn(P.vertices[i], x, P.vertices[i - 1])
                return false
            end
        end
    end
    return true
end
