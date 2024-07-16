"""
    an_element(X::Star)

Return some element of a star.

### Input

- `X` -- star

### Output

An element of the star.

### Algorithm

We apply the affine map to the result of `an_element` on the predicate.
"""
function an_element(X::Star)
    return basis(X) * an_element(predicate(X)) + center(X)
end
