"""
    vertices_list(cap::Intersection)

Return a list of vertices of a lazy intersection of two (polyhedral) sets.

### Input

- `cap` -- intersection of two (polyhedral) sets

### Output

A list containing the vertices of the lazy intersection of two sets.

### Notes

We assume that the underlying sets are polyhedral and that the intersection is
bounded.

### Algorithm

We compute the concrete intersection using `intersection` and then take the
vertices of that representation.
"""
@validate function vertices_list(cap::Intersection)
    return vertices_list(intersection(cap.X, cap.Y))
end
