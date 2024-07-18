"""
    ngens(Z::Zonotope)

Return the number of generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

An integer representing the number of generators.
"""
ngens(Z::Zonotope) = size(Z.generators, 2)
