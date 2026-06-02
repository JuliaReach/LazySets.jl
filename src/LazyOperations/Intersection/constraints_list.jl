"""
    constraints_list(cap::Intersection)

Return a list of constraints of an intersection of two (polyhedral) sets.

### Input

- `cap` -- intersection of two (polyhedral) sets

### Output

A list of constraints of the intersection.

### Notes

We assume that the underlying sets are polyhedral, i.e., offer a method
`constraints_list`.

### Algorithm

We create the polyhedron by taking the intersection of the `constraints_list`s
of the sets and remove redundant constraints.

This function ignores the boolean output from the in-place
`remove_redundant_constraints!`, which may inform the user that the constraints
are infeasible. In that case, the list of constraints at the moment when the
infeasibility was detected is returned.
"""
function constraints_list(cap::Intersection)
    constraints = [constraints_list(cap.X); constraints_list(cap.Y)]
    remove_redundant_constraints!(constraints)
    return constraints
end
