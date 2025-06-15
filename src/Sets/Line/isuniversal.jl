"""
# Extended help

    isuniversal(L::Line; [witness::Bool]=false)

### Algorithm

* If `witness` is `false`, the result is `true` if the ambient dimension is one,
and `false` otherwise.

* If `witness` is `true`, the result is `(true, [])` if the ambient dimension is
one, and `(false, v)` where ``v ∉ P`` otherwise.
"""
function isuniversal(L::Line; witness::Bool=false)
    return isuniversal(L, Val(witness))
end

# TODO implement case with witness
function isuniversal(L::Line, ::Val{false})
    return dim(L) == 1
end
