"""
    dim(P::HParallelotope)

Return the dimension of a parallelotope in constraint representation.

### Input

- `P` -- parallelotope in constraint representation

### Output

The ambient dimension of the parallelotope.
"""
function dim(P::HParallelotope)
    return size(P.directions, 1)
end
