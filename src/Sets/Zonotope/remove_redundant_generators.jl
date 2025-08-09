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
    if dim(Z) == 1  # more efficient implementation in 1D
        return _remove_redundant_generators_1d(Z)
    end
    
    G = genmat(Z)
    if isempty(G)
        return Z
    end
    G = remove_zero_columns(G)
    p = size(G, 2)

    if p <= 1
        return Zonotope(c, G)
    end

    # normalize columns 
    Gnorm = similar(G)
    @inbounds for (j, col) in enumerate(eachcol(G))
        Gnorm[:, j] = col ./ norm(col)
    end

    # sort column in ascedning order
    ord = sortperm(eachcol(Gnorm))
    merged = Vector{Vector{eltype(G)}}()

    # for each sorted column check right neighbour
    cur = copy(G[:, ord[1]])
    @inbounds for i in 2:p
        prev = ord[i-1]
        this = ord[i]

        if isapprox(view(Gnorm, :, prev), view(Gnorm, :, this))
            cur .+= G[:, this] # merge multiples
        else
            push!(merged, cur) # add column to reduced matrix 
            cur = copy(G[:, this])
        end
    end

    push!(merged, cur)
    G_new = hcat(merged...)
    
    return Zonotope(center(Z), G_new)
end

function _remove_redundant_generators_1d(Z)
    G = genmat(Z)
    g = sum(abs, G)
    return Zonotope(center(Z), hcat(g))
end
