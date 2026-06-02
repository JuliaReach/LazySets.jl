"""
    dim(cup::UnionSetArray)

Return the dimension of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

The ambient dimension of the union of a finite number of sets, or `0` if there
is no set in the array.
"""
function dim(cup::UnionSetArray)
    return length(cup.array) == 0 ? 0 : dim(cup.array[1])
end
