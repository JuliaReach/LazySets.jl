"""
    dim(P::HPoly)

Return the dimension of a polyhedron in constraint representation.

### Input

- `P`  -- polyhedron in constraint representation

### Output

The ambient dimension of the polyhedron in constraint representation.
If it has no constraints, the result is ``-1``.
"""
function dim(P::HPoly)
    return length(P.constraints) == 0 ? -1 : length(P.constraints[1].a)
end
