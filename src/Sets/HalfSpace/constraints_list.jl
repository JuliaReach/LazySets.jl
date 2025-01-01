function constraints_list(hs::HalfSpace)
    return [hs]
end

"""
    constraints_list(A::AbstractMatrix{N}, b::AbstractVector)

Convert a matrix-vector representation to a linear-constraint representation.

### Input

- `A` -- matrix
- `b` -- vector

### Output

A list of linear constraints.
"""
function constraints_list(A::AbstractMatrix, b::AbstractVector)
    m = size(A, 1)
    @assert m == length(b) "a matrix with $m rows is incompatible with a " *
                           "vector of length $(length(b))"

    return [HalfSpace(A[i, :], b[i]) for i in 1:m]
end
