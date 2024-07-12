"""
    constraints_list(L::Line2D)

Return the list of constraints of a 2D line.

### Input

- `L` -- 2D line

### Output

A list containing two half-spaces.
"""
function constraints_list(L::Line2D)
    return _constraints_list_hyperplane(L.a, L.b)
end
