export Ball1

"""
    Ball1 <: LazySet

Type that represents a ball in the 1-norm, also known as Manhattan or Taxicab
norm.

It is defined as the set

```math
\\mathcal{B}_1^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ∑_{i=1}^n |x_i| ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

### Examples

Unit ball in the 1-norm in the plane:

```jldoctest ball1_constructor
julia> B = Ball1(zeros(2), 1.)
LazySets.Ball1([0.0, 0.0], 1.)
julia> dim(B)
2
```

We evaluate the support vector in the East direction:

```jldoctest ball1_constructor
julia> σ([0.,1], B)
1-element Array{Float64,1}:
 1.0
 0.0
```
"""
struct Ball1{N<:Real} <: LazySet
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    Ball1{N}(center, radius) where N =
        radius < zero(N) ? throw(DomainError()) : new(center, radius)
end
# type-less convenience constructor
Ball1(center::Vector{N}, radius::N) where {N<:Real} = Ball1{N}(center, radius)

"""
    dim(B)

Return the dimension of a Ball1.

### Input

- `B` -- a ball in the p-norm

### Output

The ambient dimension of the ball.
"""
dim(B::Ball1) = length(B.center)

"""
    σ(d::AbstractVector{N}, B::Ball1)::AbstractVector{N} where{N<:AbstractFloat}

Return the support vector of a `Ball1` in a given direction.

### Input

- `d` -- a direction
- `B` -- a ball in the p-norm

### Output

The support vector in the given direction.

### Algorithm

The support vector of the unit ball in the ``p``-norm along direction ``d`` is:

```math
σ_{\\mathcal{B}_p^n(0, 1)}(d) = \\dfrac{\\tilde{v}}{‖\\tilde{v}‖_q},
```
where ``\\tilde{v}_i = \\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and ``tilde{v}_i = 0``
otherwise, for all ``i=1,…,n``, and ``q`` is the conjugate number of ``p``.
By the affine transformation ``x = r\\tilde{x} + c``, one obtains that
the support vector of ``\\mathcal{B}_p^n(c, r)`` is

```math
σ_{\\mathcal{B}_p^n(c, r)}(d) = \\dfrac{v}{‖v‖_q},
```
where ``v_i = c_i + r\\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and ``v_i = 0`` otherwise,
for all ``i=1,…,n``.
"""
function σ(d::AbstractVector{N}, B::Ball1)::AbstractVector{N} where{N<:AbstractFloat}
    p = B.p
    q = p/(p-1)
    v = similar(d)
    if p == one(T)
        xyz
    else
        @inbounds for (i, di) in enumerate(d)
            v[i] = di == zero(N) ? di : abs.(di).^q / di
        end
        vnorm = norm(v, p)
        svec = vnorm != zero(N) ? @.(B.center + B.radius * (v/vnorm)) : B.center
    end
    return svec
end
