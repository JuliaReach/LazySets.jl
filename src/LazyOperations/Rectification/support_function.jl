"""
    ρ(d::AbstractVector, R::Rectification)

Evaluate the support function of a rectification in a given direction.

### Input

- `d` -- direction
- `R` -- rectification

### Output

The support value of the rectification in the given direction.

### Algorithm

We use different procedures for different types of input sets.
If the wrapped set has a suitable structure for which we can efficiently compute
the support vector, we fall back to the evaluation of the support function by
means of the support vector.
Otherwise we compute the union of projections to obtain a precise result (see
[`to_union_of_projections`](@ref)), and then compute the support function for
this union.
(The union is cached internally, so subsequent queries are more efficient.)
"""
@validate function ρ(d::AbstractVector, R::Rectification)
    if R.cache.use_support_vector
        return dot(d, σ(d, R))
    end

    _compute_exact_representation!(R)
    return ρ(d, R.cache.set)
end
