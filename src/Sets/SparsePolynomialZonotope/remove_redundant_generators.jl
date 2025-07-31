"""
# Extended help

    remove_redundant_generators(S::SparsePolynomialZonotope)

## Notes

The result uses dense arrays irrespective of the array type of `S`.

### Algorithm

Let `G` be the dependent generator matrix, `E` the exponent matrix, and `GI` the
independent generator matrix of `S`. We perform the following simplifications:

- Remove zero columns in `G` and the corresponding columns in `E`.
- Remove Zero columns in `GI`.
- For zero columns in `E`, add the corresponding column in `G` to the center.
- Group repeated columns in `E` together by summing the corresponding columns in
  `G`.
"""
function remove_redundant_generators(S::SparsePolynomialZonotope)
    c, G, E = _remove_redundant_generators_polyzono(center(S), genmat_dep(S),
                                                    expmat(S))
    GI = remove_zero_columns(genmat_indep(S))
    return SparsePolynomialZonotope(c, G, GI, E, indexvector(S))
end
