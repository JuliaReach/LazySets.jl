"""
    generators(P::HParallelotope)

Return an iterator over the generators of a parallelotope in constraint
representation.

### Input

- `P` -- parallelotope in constraint representation

### Output

An iterator over the generators of `P`.
"""
function generators(P::HParallelotope)
    return generators_fallback(P)
end
