"""
    dim(ia::IntersectionArray)

Return the dimension of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of sets

### Output

The ambient dimension of the intersection of a finite number of sets, or `0` if
there is no set in the array.
"""
function dim(ia::IntersectionArray)
    return length(ia.array) == 0 ? 0 : dim(ia.array[1])
end
