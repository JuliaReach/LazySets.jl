"""
    addconstraint!(P::HPoly, constraint::HalfSpace)

Add a linear constraint to a polyhedron in constraint representation.

### Input

- `P`          -- polyhedron in constraint representation
- `constraint` -- linear constraint to add

### Notes

It is left to the user to guarantee that the dimension of all linear constraints
is the same.
"""
function addconstraint!(P::HPoly, constraint::HalfSpace)
    push!(P.constraints, constraint)
    return nothing
end
