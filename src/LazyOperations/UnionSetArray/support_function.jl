"""
    ρ(d::AbstractVector, cup::UnionSetArray)

Evaluate the support function of the union of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `cup` -- union of a finite number of sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function of the union of a finite number of sets ``X₁, X₂, ...``
can be obtained as the maximum of ``ρ(d, X₂), ρ(d, X₂), ...``.
"""
@validate function ρ(d::AbstractVector, cup::UnionSetArray)
    return maximum(Xi -> ρ(d, Xi), array(cup))
end
