"""
    dim(C::Complement)

Return the dimension of the complement of a set.

### Input

- `C` -- complement of a set

### Output

The ambient dimension of the complement of a set.
"""
function dim(C::Complement)
    return dim(C.X)
end
