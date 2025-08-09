"""
# Extended help

    remove_redundant_generators(Z::Zonotope)

### Algorithm

This method normalizes each generator to have unit norm, sorts the normalized 
generators lexicographicaly, and then traverses them once to merge consecutive 
generators that are approximately collinear by summing them into a single generator.

For one-dimensional zonotopes we use a more efficient implementation where we
just take the absolute sum of all generators.
"""
function remove_redundant_generators(Z::Zonotope)
    G = genmat(Z)
    Gnew = remove_redundant_generators(G)

    if Gnew === G
        return Z  # return the original zonotope if no generator was removed
    end
    return Zonotope(center(Z), Gnew)
end
