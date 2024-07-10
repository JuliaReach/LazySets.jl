"""
    translate!(L::Line, v::AbstractVector)

Translate (i.e., shift) a line by a given vector, in-place.

### Input

- `L` -- line
- `v` -- translation vector

### Output

The line `L` translated by `v`.
"""
function translate!(L::Line, v::AbstractVector)
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"

    L.p .+= v
    return L
end
