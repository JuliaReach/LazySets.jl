"""
    in(x::AbstractVector, lm::LinearMap)

Check whether a given point is contained in a linear map.

### Input

- `x`  -- point/vector
- `lm` -- linear map

### Output

`true` iff ``x ∈ lm``.

### Algorithm

For well-conditioned square or tall matrices, we compute ``y = M \\ x`` and
check whether ``y ∈ S``. For tall matrices, we also check ``M⋅y = x``.
For wide or ill-conditioned matrices, we use a concrete representation of the
linear map and check membership in that set.

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
    M = matrix(lm)
    if size(M, 1) < size(M, 2) || !iswellconditioned(M)
        # wide or ill-conditioned matrix; use concrete set representation
        return x ∈ linear_map(M, set(lm))
    end
    y = M \ x
    return (issquare(M) || _isapprox(M * y, x)) && y ∈ set(lm)
end
