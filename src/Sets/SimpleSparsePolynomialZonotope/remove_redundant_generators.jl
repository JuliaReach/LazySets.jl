"""
    remove_redundant_generators(S::SimpleSparsePolynomialZonotope)

Remove redundant generators from a simple sparse polynomial zonotope.

### Input

- `S` -- simple sparse polynomial zonotope

### Output

A new simple sparse polynomial zonotope such that redundant generators have been
removed.

## Notes

The result uses dense arrays irrespective of the array type of `S`.

### Algorithm

Let `G` be the generator matrix and `E` the exponent matrix of `S`. The
following simplifications are performed:

- Zero columns in `G` and the corresponding columns in `E` are removed.
- For zero columns in `E`, the corresponding column in `G` is summed to the
  center.
- Repeated columns in `E` are grouped together by summing the corresponding
  columns in `G`.
"""
function remove_redundant_generators(S::SimpleSparsePolynomialZonotope)
    c, G, E = _remove_redundant_generators_polyzono(center(S), genmat(S), expmat(S))

    return SimpleSparsePolynomialZonotope(c, G, E)
end

function _remove_redundant_generators_polyzono(c, G, E)
    Gnew = Matrix{eltype(G)}(undef, size(G, 1), 0)
    Enew = Matrix{eltype(E)}(undef, size(E, 1), 0)
    cnew = copy(c)

    visited_exps = Dict{Vector{Int},Int}()
    @inbounds for (gi, ei) in zip(eachcol(G), eachcol(E))
        all(isapproxzero, gi) && continue
        if iszero(ei)
            cnew += gi
        elseif haskey(visited_exps, ei) # repeated exponent
            idx = visited_exps[ei]
            Gnew[:, idx] += gi
        else
            Gnew = hcat(Gnew, gi)
            Enew = hcat(Enew, ei)
            visited_exps[ei] = size(Enew, 2)
        end
    end

    return cnew, Gnew, Enew
end
