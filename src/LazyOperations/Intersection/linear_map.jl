"""
    linear_map(M::AbstractMatrix, cap::Intersection)

Return the concrete linear map of an intersection of two sets.

### Input

- `M`   -- matrix
- `cap` -- intersection of two sets

### Output

The set obtained by applying the given linear map to the intersection.

### Algorithm

This method computes the concrete intersection.
"""
@validate function linear_map(M::AbstractMatrix, cap::Intersection)
    return linear_map(M, intersection(cap.X, cap.Y))
end
