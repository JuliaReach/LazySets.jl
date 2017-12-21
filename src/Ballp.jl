import Base: ∈, ⊆

export Ballp

"""
    Ballp{N<:AbstractFloat} <: AbstractPointSymmetric{N}

Type that represents a ball in the p-norm, for ``1 ≤ p ≤ ∞``.

It is defined as the set

```math
\\mathcal{B}_p^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_p ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_p`` for ``1 ≤ p ≤ ∞`` denotes the vector ``p``-norm, defined as
``‖ x ‖_p = \\left( \\sum\\limits_{i=1}^n |x_i|^p \\right)^{1/p}`` for any
``x ∈ \\mathbb{R}^n``.

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
LazySets.Ballp{Float64}(1.5, [0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
julia> dim(B)
5
```

We evaluate the support vector in direction ``[1,2,…,5]``:

```jldoctest ballp_constructor
julia> σ(1.:5, B)
5-element Array{Float64,1}:
 0.013516
 0.054064
 0.121644
 0.216256
 0.3379
```
"""
struct Ballp{N<:AbstractFloat} <: AbstractPointSymmetric{N}
    p::N
    center::Vector{N}
    radius::N

    function Ballp{N}(p, center, radius) where N
        if radius < zero(N)
            throw(DomainError())
        end
        if p == Inf
            return BallInf(center, radius)
        elseif p == 2
            return Ball2(center, radius)
        elseif p == 1
            return Ball1(center, radius)
        elseif 1 < p && p < Inf
            new(p, center, radius)
        else
            throw(DomainError())
        end
    end
end
# type-less convenience constructor
Ballp(p::N, center::Vector{N}, radius::N) where {N<:AbstractFloat} =
    Ballp{N}(p, center, radius)


# --- AbstractPointSymmetric interface functions ---


"""
    center(B::Ballp{N})::Vector{N} where {N<:AbstractFloat}

Return the center of a ball in the p-norm.

### Input

- `B` -- ball in the p-norm

### Output

The center of the ball in the p-norm.
"""
function center(B::Ballp{N})::Vector{N} where {N<:AbstractFloat}
    return B.center
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, B::Ballp)::AbstractVector{N} where {N<:AbstractFloat}

Return the support vector of a `Ballp` in a given direction.

### Input

- `d` -- direction
- `B` -- ball in the p-norm

### Output

The support vector in the given direction.
If the direction has norm zero, the center of the ball is returned.

### Algorithm

The support vector of the unit ball in the ``p``-norm along direction ``d`` is:

```math
σ_{\\mathcal{B}_p^n(0, 1)}(d) = \\dfrac{\\tilde{v}}{‖\\tilde{v}‖_q},
```
where ``\\tilde{v}_i = \\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and
``tilde{v}_i = 0`` otherwise, for all ``i=1,…,n``, and ``q`` is the conjugate
number of ``p``.
By the affine transformation ``x = r\\tilde{x} + c``, one obtains that
the support vector of ``\\mathcal{B}_p^n(c, r)`` is

```math
σ_{\\mathcal{B}_p^n(c, r)}(d) = \\dfrac{v}{‖v‖_q},
```
where ``v_i = c_i + r\\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and ``v_i = 0``
otherwise, for all ``i = 1, …, n``.
"""
function σ(d::AbstractVector{N},
           B::Ballp)::AbstractVector{N} where {N<:AbstractFloat}
    p = B.p
    q = p/(p-1)
    v = similar(d)
    @inbounds for (i, di) in enumerate(d)
        v[i] = di == zero(N) ? di : abs.(di).^q / di
    end
    vnorm = norm(v, p)
    svec = vnorm != zero(N) ? @.(B.center + B.radius * (v/vnorm)) : B.center
    return svec
end

"""
    ∈(x::AbstractVector{N}, B::Ballp{N})::Bool where {N<:AbstractFloat}

Check whether a given point is contained in a ball in the p-norm.

### Input

- `x` -- point/vector
- `B` -- ball in the p-norm

### Output

`true` iff ``x ∈ B``.

### Notes

This implementation is worst-case optimized, i.e., it is optimistic and first
computes (see below) the whole sum before comparing to the radius.
In applications where the point is typically far away from the ball, a fail-fast
implementation with interleaved comparisons could be more efficient.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the p-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``\\left( ∑_{i=1}^n |c_i - x_i|^p \\right)^{1/p} ≤ r``.

### Examples

```jldoctest
julia> B = Ballp(1.5, [1., 1.], 1.)
LazySets.Ballp{Float64}(1.5, [1.0, 1.0], 1.0)
julia> ∈([.5, -.5], B)
false
julia> ∈([.5, 1.5], B)
true
```
"""
function ∈(x::AbstractVector{N}, B::Ballp{N})::Bool where {N<:AbstractFloat}
    @assert length(x) == dim(B)
    sum = zero(N)
    for i in eachindex(x)
        sum += abs(B.center[i] - x[i])^B.p
    end
    return sum^(1./B.p) <= B.radius
end

"""
    ⊆(B::Ballp, S::AbstractSingleton)::Bool

Check whether a ball in the p-norm is contained in a set with a single value.

### Input

- `B` -- inner ball in the p-norm
- `S` -- outer set with a single value

### Output

`true` iff ``B ⊆ S``.
"""
function ⊆(B::Ballp, S::AbstractSingleton)::Bool
    return B.center == element(S) && B.radius == 0
end

"""
    ⊆(B::Ballp, H::AbstractHyperrectangle)::Bool

Check whether a ball in the p-norm is contained in a hyperrectangle.

### Input

- `B` -- inner ball in the p-norm
- `H` -- outer hyperrectangle

### Output

`true` iff ``B ⊆ H``.

### Algorithm

This implementation computes the interval hull of the ball and then checks
containment in the hyperrectangle.
"""
function ⊆(B::Ballp, H::AbstractHyperrectangle)::Bool
    return ⊆(Approximations.box_approximation(B), H)
end
