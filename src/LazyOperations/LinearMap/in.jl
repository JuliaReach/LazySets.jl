"""
    in(x::AbstractVector, lm::LinearMap)

Check whether a given point is contained in a linear map.

### Input

- `x`  -- point/vector
- `lm` -- linear map

### Output

`true` iff ``x ∈ lm``.

### Algorithm

Note that ``x ∈ M⋅S`` iff ``M^{-1}⋅x ∈ S``.
This implementation does not explicitly invert the matrix: instead of
``M^{-1}⋅x`` it computes ``M \\ x``.
Hence it also works for non-square matrices.

### Examples

```jldoctest
julia> lm = LinearMap([2.0 0.0; 0.0 1.0], BallInf([1., 1.], 1.));

julia> [5.0, 1.0] ∈ lm
false
julia> [3.0, 1.0] ∈ lm
true
```

An example with non-square matrix:
```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> [0.5, 0.5] ∈ M*B
true
```
"""
@validate function in(x::AbstractVector, lm::LinearMap)
    if !iswellconditioned(matrix(lm))
        # ill-conditioned matrix; use concrete set representation
        return x ∈ linear_map(matrix(lm), set(lm))
    end
    return matrix(lm) \ x ∈ set(lm)
end
