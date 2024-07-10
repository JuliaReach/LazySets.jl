"""
    isuniversal(L::Line; [witness::Bool]=false)

Check whether a line is universal.

### Input

- `P`       -- line
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` is `false`: `true` if the ambient dimension is one, `false`
otherwise.

* If `witness` is `true`: `(true, [])` if the ambient dimension is one,
`(false, v)` where ``v ∉ P`` otherwise.
"""
isuniversal(L::Line; witness::Bool=false) = isuniversal(L, Val(witness))

# TODO implement case with witness
isuniversal(L::Line, ::Val{false}) = dim(L) == 1
