"""
# Extended help

    remove_redundant_generators(Z::Zonotope)

### Algorithm

For each generator ``g_j`` that has not been checked yet, we find all other
generators that are linearly dependent with ``g_j``.
Then we combine those generators into a single generator.

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
