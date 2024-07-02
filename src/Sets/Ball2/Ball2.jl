"""
    Ball2{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractBallp{N}

Type that represents a ball in the 2-norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Notes

Mathematically, a ball in the 2-norm is defined as the set

```math
\\mathcal{B}_2^n(c, r) = \\{ x ∈ ℝ^n : ‖ x - c ‖_2 ≤ r \\},
```
where ``c ∈ ℝ^n`` is its center and ``r ∈ ℝ_+`` its radius.
Here ``‖ ⋅ ‖_2`` denotes the Euclidean norm (also known as 2-norm), defined as
``‖ x ‖_2 = \\left( ∑\\limits_{i=1}^n |x_i|^2 \\right)^{1/2}`` for any
``x ∈ ℝ^n``.

### Examples

Create a five-dimensional ball `B` in the 2-norm centered at the origin with
radius 0.5:

```jldoctest ball2_label
julia> B = Ball2(zeros(5), 0.5)
Ball2{Float64, Vector{Float64}}([0.0, 0.0, 0.0, 0.0, 0.0], 0.5)

julia> dim(B)
5
```

Evaluate `B`'s support vector in the direction ``[1,2,3,4,5]``:

```jldoctest ball2_label
julia> σ([1.0, 2, 3, 4, 5], B)
5-element Vector{Float64}:
 0.06741998624632421
 0.13483997249264842
 0.20225995873897262
 0.26967994498529685
 0.3370999312316211
```
"""
struct Ball2{N<:AbstractFloat,VN<:AbstractVector{N}} <: AbstractBallp{N}
    center::VN
    radius::N

    # default constructor with domain constraint for radius
    function Ball2(center::VN, radius::N) where {N<:AbstractFloat,VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must not be negative"
        return new{N,VN}(center, radius)
    end
end
