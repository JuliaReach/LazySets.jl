"""
    tosimplehrep(constraints::AbstractVector{<:HalfSpace}; [n]::Int=0)

Return the simple H-representation ``Ax ≤ b`` from a list of linear constraints.

### Input

- `constraints` -- a list of linear constraints
- `n`           -- (optional; default: `0`) dimension of the constraints

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.

### Notes

The parameter `n` can be used to create a matrix with no constraints but a
non-zero dimension. If `n` ≤ 0, this method takes the dimension of the first
constraint.
"""
function tosimplehrep(constraints::AbstractVector{<:HalfSpace}; n::Int=0)
    N = eltype(eltype(constraints))
    m = length(constraints)
    if m == 0
        A = Matrix{N}(undef, 0, n)
        b = Vector{N}(undef, 0)
        return (A, b)
    end
    if n <= 0
        n = dim(first(constraints))
    end
    A = zeros(N, m, n)
    b = zeros(N, m)
    @inbounds begin
        for (i, Pi) in enumerate(constraints)
            A[i, :] = Pi.a
            b[i] = Pi.b
        end
    end
    return (A, b)
end
