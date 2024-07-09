"""
    translate!(S::Singleton, v::AbstractVector)

Translate (i.e., shift) a singleton by a given vector, in-place.

### Input

- `S` -- singleton
- `v` -- translation vector

### Output

The singleton `S` translated by `v`.

### Algorithm

We add the vector to the point in the singleton.

### Notes

See also [`translate(::Singleton, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(S::Singleton, v::AbstractVector)
    @assert length(v) == dim(S) "cannot translate a $(dim(S))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    S.element .+= v
    return S
end
