"""
    constraints_list(ia::IntersectionArray)

Return the list of constraints of an intersection of a finite number of
(polyhedral) sets.

### Input

- `ia` -- intersection of a finite number of (polyhedral) sets

### Output

The list of constraints of the intersection.

### Notes

We assume that the underlying sets are polyhedral, i.e., offer a method
`constraints_list`.

### Algorithm

We create the polyhedron from the `constraints_list`s of the sets and remove
redundant constraints.
"""
function constraints_list(ia::IntersectionArray)
    N = eltype(ia)
    constraints = Vector{HalfSpace{N,Vector{N}}}() # TODO: use vector type of ia
    for X in ia
        clist_X = _normal_Vector(X)
        if clist_X isa HalfSpace
            push!(constraints, clist_X)
        else
            append!(constraints, clist_X)
        end
    end
    remove_redundant_constraints!(constraints)
    return constraints
end
