"""
# Extended help

    isuniversal(L::Line; [witness::Bool]=false)

### Algorithm

* If `witness` is `false`, the result is `true` if the ambient dimension is one,
and `false` otherwise.

* If `witness` is `true`, the result is `(true, [])` if the ambient dimension is
one, and `(false, v)` where ``v ∉ P`` otherwise.
"""
isuniversal(L::Line; witness::Bool=false) = isuniversal(L, Val(witness))

# TODO implement case with witness
isuniversal(L::Line, ::Val{false}) = dim(L) == 1
