"""
    Ballp{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractBallp{N}

Type that represents a ball in the p-norm, for ``1 ≤ p ≤ ∞``.

It is defined as the set

```math
\\mathcal{B}_p^n(c, r) = \\{ x ∈ ℝ^n : ‖ x - c ‖_p ≤ r \\},
```
where ``c ∈ ℝ^n`` is its center and ``r ∈ ℝ_+`` its radius.
Here ``‖ ⋅ ‖_p`` for ``1 ≤ p ≤ ∞`` denotes the vector ``p``-norm, defined as
``‖ x ‖_p = \\left( ∑\\limits_{i=1}^n |x_i|^p \\right)^{1/p}`` for any
``x ∈ ℝ^n``.

### Fields

- `p`      -- norm as a real scalar
- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

## Notes

The special cases ``p=1``, ``p=2`` and ``p=∞`` fall back to the specialized
types `Ball1`, `Ball2` and `BallInf`, respectively.

### Examples

A five-dimensional ball in the ``p=3/2`` norm centered at the origin of radius
0.5:

```jldoctest ballp_constructor
julia> B = Ballp(3/2, zeros(5), 0.5)
Ballp{Float64, Vector{Float64}}(1.5, [0.0, 0.0, 0.0, 0.0, 0.0], 0.5)

julia> dim(B)
5
```

We evaluate the support vector in direction ``[1,2,…,5]``:

```jldoctest ballp_constructor
julia> σ([1.0, 2, 3, 4, 5], B)
5-element Vector{Float64}:
 0.013516004434607558
 0.05406401773843023
 0.12164403991146802
 0.21625607095372093
 0.33790011086518895
```
"""
struct Ballp{N<:AbstractFloat,VN<:AbstractVector{N}} <: AbstractBallp{N}
    p::N
    center::VN
    radius::N

    # default constructor with domain constraint for radius and p
    function Ballp(p::N, center::VN, radius::N) where {N,VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must be nonnegative but is $radius"
        @assert isfinite(radius) "the radius must be finite but is $radius"
        @assert p >= one(N) "p must not be less than 1"
        if p == N(Inf)
            require(@__MODULE__, :LazySets; fun_name="Ballp")

            return BallInf(center, radius)
        elseif p == N(2)
            require(@__MODULE__, :LazySets; fun_name="Ballp")

            return Ball2(center, radius)
        elseif isone(p)
            require(@__MODULE__, :LazySets; fun_name="Ballp")

            return Ball1(center, radius)
        else
            return new{N,VN}(p, center, radius)
        end
    end
end
