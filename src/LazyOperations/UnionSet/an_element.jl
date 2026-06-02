"""
    an_element(cup::UnionSet)

Return some element of the union of two sets.

### Input

- `cup` -- union of two sets

### Output

An element in the union of two sets.

### Algorithm

We use `an_element` on the first non-empty wrapped set.
"""
function an_element(cup::UnionSet)
    if isempty(cup.X)
        return an_element(cup.Y)
    end
    return an_element(cup.X)
end
