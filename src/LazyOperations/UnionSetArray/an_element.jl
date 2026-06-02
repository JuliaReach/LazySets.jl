"""
    an_element(cup::UnionSetArray)

Return some element of the union of a finite number of sets.

### Input

- `cup` -- union of a finite number of sets

### Output

An element in the union of a finite number of sets.

### Algorithm

We use `an_element` on the first non-empty wrapped set.
"""
function an_element(cup::UnionSetArray)
    for Xi in cup
        if !isempty(Xi)
            return an_element(Xi)
        end
    end
    return throw(ArgumentError("an empty set does not have any element"))
end
