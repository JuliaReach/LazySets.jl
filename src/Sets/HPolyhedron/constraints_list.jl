"""
    constraints_list(P::HPoly)

Return the list of constraints defining a polyhedron in constraint
representation.

### Input

- `P` -- polyhedron in constraint representation

### Output

The list of constraints of the polyhedron.
"""
function constraints_list(P::HPoly)
    return P.constraints
end
