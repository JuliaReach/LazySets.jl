"""
# Extended help

    remove_redundant_generators(Z::Zonotope)

### Algorithm

This method normalizes each generator to have unit norm, sorts the normalized 
generators lexicographicaly, and then traverses them once to merge consecutive 
generators that are approximately collinear by summing them into a single generator.

For `Rational` inputs, we use a less efficient implementation: For each generator ``g_j`` that has
not been checked yet, we find all other generators that are linearly dependent with ``g_j``.
Then we combine those generators into a single generator.

For one-dimensional zonotopes we use a more efficient implementation where we
just take the absolute sum of all generators.

### Notes

The implementation may change the order of the non-removed generators.
"""
function remove_redundant_generators(Z::Zonotope)
    G = genmat(Z)
    Gnew = remove_redundant_generators(G)

    if Gnew === G
        return Z  # return the original zonotope if no generator was removed
    end
    return Zonotope(center(Z), Gnew)
end
