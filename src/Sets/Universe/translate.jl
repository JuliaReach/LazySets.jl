"""
    translate(U::Universe, v::AbstractVector)

Translate (i.e., shift) a universe by a given vector.

### Input

- `U` -- universe
- `v` -- translation vector

### Output

The universe.
"""
function translate(U::Universe, v::AbstractVector)
    return translate!(U, v)  # no need to copy
end

"""
    translate!(U::Universe, v::AbstractVector)

Translate (i.e., shift) a universe by a given vector, in-place.

### Input

- `U` -- universe
- `v` -- translation vector

### Output

The universe.
"""
function translate!(U::Universe, v::AbstractVector)
    @assert length(v) == dim(U) "cannot translate a $(dim(U))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return U
end
