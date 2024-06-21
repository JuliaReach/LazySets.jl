"""
    translate(∅::EmptySet, v::AbstractVector)

Translate (i.e., shift) an empty set by a given vector.

### Input

- `∅` -- empty set
- `v` -- translation vector

### Output

The empty set.
"""
function translate(∅::EmptySet, v::AbstractVector)
    return translate!(∅, v)  # no need to copy
end

function translate!(∅::EmptySet, v::AbstractVector)
    @assert length(v) == dim(∅) "cannot translate a $(dim(∅))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return ∅
end
