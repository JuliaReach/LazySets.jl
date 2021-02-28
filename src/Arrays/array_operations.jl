export same_sign

"""
    same_sign(A::AbstractArray{N}; [optimistic]::Bool=false) where {N}

Check whether all elements of the given array have the same sign.

### Input

- `A`          -- array
- `optimistic` -- (optional; default: `false`) flag for expressing that the
                  expected result is `true`

### Output

`true` if and only if all elements in `M` have the same sign.

### Algorithm

If `optimistic` is `false`, we check the sign of the first element and compare
to the sign of all elements.

If `optimistic` is `true`, we compare the absolute element sum with the sum of
the absolute of the elements; this is faster if the result is `true` because
there is no branching.

```math
    |\\sum_i A_i| = \\sum_i |A_i|
```
"""
function same_sign(A::AbstractArray{N}; optimistic::Bool=false) where {N}
    if optimistic
        return sum(abs, A) == abs(sum(A))
    else
        if isempty(A)
            return true
        end
        @inbounds if first(A) >= zero(N)
            return all(e -> e >= zero(N), A)
        else
            return all(e -> e <= zero(N), A)
        end
    end
end
