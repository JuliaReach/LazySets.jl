"""
    BallInf{N, VN<:AbstractVector{N}} <: AbstractHyperrectangle{N}

Type that represents a ball in the infinity norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Notes

Mathematically, a ball in the infinity norm is defined as the set

```math
\\mathcal{B}_∞^n(c, r) = \\{ x ∈ ℝ^n : ‖ x - c ‖_∞ ≤ r \\},
```
where ``c ∈ ℝ^n`` is its center and ``r ∈ ℝ_+`` its radius.
Here ``‖ ⋅ ‖_∞`` denotes the infinity norm, defined as
``‖ x ‖_∞ = \\max\\limits_{i=1,…,n} | x_i |`` for any
``x ∈ ℝ^n``.

### Examples

Construct the two-dimensional unit ball and compute its support function along
the positive ``x=y`` direction:

```jldoctest
julia> B = BallInf(zeros(2), 1.0)
BallInf{Float64, Vector{Float64}}([0.0, 0.0], 1.0)

julia> dim(B)
2

julia> ρ([1.0, 1.0], B)
2.0
```
"""
struct BallInf{N,VN<:AbstractVector{N}} <: AbstractHyperrectangle{N}
    center::VN
    radius::N

    # default constructor with domain constraint for radius
    function BallInf(center::VN, radius::N) where {N,VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must be nonnegative but is $radius"
        return new{N,VN}(center, radius)
    end
end
