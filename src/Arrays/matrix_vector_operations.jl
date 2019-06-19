export inner

"""
    inner(x::AbstractVector{N}, A::AbstractMatrix{N}, y::AbstractVector{N}
         ) where {N}

Compute the inner product ``xᵀ A y``.

### Input

- `x` -- vector on the left
- `A` -- matrix
- `y` -- vector on the right

### Output

The (scalar) result of the multiplication.
"""
function inner(x::AbstractVector{N}, A::AbstractMatrix{N}, y::AbstractVector{N}
              ) where {N}
    dot(x, A * y)
end
