"""
    directions(P::HParallelotope)

Return the directions matrix of a parallelotope in constraint representation.

### Input

- `P` -- parallelotope in constraint representation

### Output

A matrix where each row represents a direction of the parallelotope.
The negated directions `-D_i` are implicit (see [`HParallelotope`](@ref) for
details).
"""
function directions(P::HParallelotope)
    return P.directions
end
