"""
    translate(P::VPolytope, v::AbstractVector)

Translate (i.e., shift) a polytope in vertex representation by a given vector.

### Input

- `P` -- polytope in vertex representation
- `v` -- translation vector

### Output

A translated polytope in vertex representation.

### Notes

See [`translate!(::VPolytope, ::AbstractVector)`](@ref) for the in-place version.
"""
function translate(P::VPolytope, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

"""
    translate!(P::VPolytope, v::AbstractVector)

Translate (i.e., shift) a polytope in vertex representation by a given vector,
in-place.

### Input

- `P` -- polytope in vertex representation
- `v` -- translation vector

### Output

The polytope `P` translated by `v`.

### Notes

See [`translate(::VPolytope, ::AbstractVector)`](@ref) for the out-of-place
version.

### Algorithm

We add the vector to each vertex of the polytope.
"""
function translate!(P::VPolytope, v::AbstractVector)
    if isempty(P.vertices)
        return P
    end

    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    for x in P.vertices
        x .+= v
    end
    return P
end
