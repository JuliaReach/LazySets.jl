"""
    remove_redundant_generators(Z::Zonotope{N}) where {N}

Remove all redundant (pairwise linearly dependent) generators of a zonotope.

### Input

- `Z` -- zonotope

### Output

A new zonotope with fewer generators, or the same zonotope if no generator could
be removed.

### Algorithm

For each generator ``g_j`` that has not been checked yet, we find all other
generators that are linearly dependent with ``g_j``.
Then we combine those generators into a single generator.

For one-dimensional zonotopes we use a more efficient implementation where we
just take the absolute sum of all generators.
"""
function remove_redundant_generators(Z::Zonotope{N}) where {N}
    if dim(Z) == 1  # more efficient implementation in 1D
        return _remove_redundant_generators_1d(Z)
    end

    G = genmat(Z)
    G = remove_zero_columns(G)
    p = size(G, 2)
    removed_zero_generators = p < ngens(Z)
    deleted = false
    done = falses(p)
    G_new = vector_type(typeof(G))[]  # list of new column vectors
    @inbounds for j1 in 1:p
        if done[j1]  # skip if the generator was already removed
            continue
        end
        # "done[j1] = true" not needed because we will never look at it again
        gj1 = G[:, j1]
        for j2 in (j1 + 1):p  # look at all generators to the right
            if done[j2]  # skip if the generator was already removed
                continue
            end
            gj2 = G[:, j2]
            answer, factor = ismultiple(gj1, gj2)
            if answer
                # column j2 is a multiple of column j1
                if factor > zero(N)
                    gj1 += gj2
                else
                    gj1 -= gj2
                end
                done[j2] = true
                deleted = true
            end
        end
        push!(G_new, gj1)
    end

    if deleted
        G_new = reduce(hcat, G_new)  # convert list of column vectors to matrix
        return Zonotope(center(Z), G_new)
    elseif removed_zero_generators
        return Zonotope(center(Z), G)
    end
    return Z  # return the original zonotope if no generator was removed
end

function _remove_redundant_generators_1d(Z)
    G = genmat(Z)
    g = sum(abs, G)
    return Zonotope(center(Z), hcat(g))
end
