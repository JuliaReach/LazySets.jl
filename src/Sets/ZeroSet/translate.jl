"""
    translate(Z::ZeroSet, v::AbstractVector)

Translate (i.e., shift) a zero set by a given vector.

### Input

- `Z` -- zero set
- `v` -- translation vector

### Output

A singleton containing the vector `v`.
"""
function translate(Z::ZeroSet, v::AbstractVector)
    require(@__MODULE__, :LazySets; fun_name="translate")

    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"

    return Singleton(v)
end
