import Base.∈

export Ball2,
       is_intersection_empty

"""
    Ball2{N<:AbstractFloat} <: AbstractPointSymmetric{N}

Type that represents a ball in the 2-norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Notes

Mathematically, a ball in the 2-norm is defined as the set

```math
\\mathcal{B}_2^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_2 ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_2`` denotes the Euclidean norm (also known as 2-norm), defined as
``‖ x ‖_2 = \\left( \\sum\\limits_{i=1}^n |x_i|^2 \\right)^{1/2}`` for any
``x ∈ \\mathbb{R}^n``.

### Examples

Create a five-dimensional ball `B` in the 2-norm centered at the origin with
radius 0.5:

```jldoctest ball2_label
julia> B = Ball2(zeros(5), 0.5)
LazySets.Ball2{Float64}([0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
julia> dim(B)
5
```

Evaluate `B`'s support vector in the direction ``[1,2,3,4,5]``:

```jldoctest ball2_label
julia> σ([1.,2.,3.,4.,5.], B)
5-element Array{Float64,1}:
 0.06742
 0.13484
 0.20226
 0.26968
 0.3371
```
"""
struct Ball2{N<:AbstractFloat} <: AbstractPointSymmetric{N}
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    Ball2{N}(center, radius) where N =
        radius < zero(N) ? throw(DomainError()) : new(center, radius)
end
# type-less convenience constructor
Ball2(center::Vector{N}, radius::N) where {N<:AbstractFloat} =
    Ball2{N}(center, radius)


# --- AbstractPointSymmetric interface functions ---


"""
    center(B::Ball2{N})::Vector{N} where {N<:AbstractFloat}

Return the center of a ball in the 2-norm.

### Input

- `B` -- ball in the 2-norm

### Output

The center of the ball in the 2-norm.
"""
function center(B::Ball2{N})::Vector{N} where {N<:AbstractFloat}
    return B.center
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, B::Ball2)::AbstractVector{<:AbstractFloat} where {N<:AbstractFloat}

Return the support vector of a 2-norm ball in a given direction.

### Input

- `d` -- direction
- `B` -- ball in the 2-norm

### Output

The support vector in the given direction.
If the direction has norm zero, the origin is returned.

### Notes

Let ``c`` and ``r`` be the center and radius of a ball ``B`` in the 2-norm,
respectively.
For nonzero direction ``d`` we have

```math
σ_B(d) = c + r \\frac{d}{‖d‖_2}.
```

This function requires computing the 2-norm of the input direction, which is
performed in the given precision of the numeric datatype of both the direction
and the set.
Exact inputs are not supported.
"""
function σ(d::AbstractVector{N},
           B::Ball2)::AbstractVector{<:AbstractFloat} where {N<:AbstractFloat}
    dnorm = norm(d, 2)
    if dnorm <= zero(N)
        return zeros(eltype(d), length(d))
    else
        return @. B.center + d * (B.radius / dnorm)
    end
end

"""
    ∈(x::AbstractVector{N}, B::Ball2{N})::Bool where {N<:AbstractFloat}

Check whether a given point is contained in a ball in the 2-norm.

### Input

- `x` -- point/vector
- `B` -- ball in the 2-norm

### Output

`true` iff ``x ∈ B``.

### Notes

This implementation is worst-case optimized, i.e., it is optimistic and first
computes (see below) the whole sum before comparing to the radius.
In applications where the point is typically far away from the ball, a fail-fast
implementation with interleaved comparisons could be more efficient.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the 2-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``\\left( ∑_{i=1}^n |c_i - x_i|^2 \\right)^{1/2} ≤ r``.

### Examples

```jldoctest
julia> B = Ball2([1., 1.], sqrt(0.5))
LazySets.Ball2{Float64}([1.0, 1.0], 0.7071067811865476)
julia> ∈([.5, 1.6], B)
false
julia> ∈([.5, 1.5], B)
true
```
"""
function ∈(x::AbstractVector{N}, B::Ball2{N})::Bool where {N<:AbstractFloat}
    @assert length(x) == dim(B)
    sum = zero(N)
    for i in eachindex(x)
        sum += (B.center[i] - x[i])^2
    end
    return sqrt(sum) <= B.radius
end


"""
    is_intersection_empty(B1::Ball2{N},
                          B2::Ball2{N},
                          witness::Bool=false
                         )::Union{Bool, Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether two balls in the 2-norm do not intersect, and otherwise optionally
compute a witness.

### Input

- `B1` -- first ball in the 2-norm
- `B2` -- second ball in the 2-norm
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``B1 ∩ B2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``B1 ∩ B2 = ∅``
  * `(false, v)` iff ``B1 ∩ B2 ≠ ∅`` and ``v ∈ B1 ∩ B2``

### Algorithm

``B1 ∩ B2 = ∅`` iff ``‖ c_2 - c_1 ‖_2 > r_1 + r_2``.

A witness is computed depending on the smaller/bigger ball (to break ties,
choose `B1` for the smaller ball) as follows.
- If the smaller ball's center is contained in the bigger ball, we return it.
- Otherwise start in the smaller ball's center and move toward the other center
  until hitting the smaller ball's border.
  In other words, the witness is the point in the smaller ball that is closest
  to the center of the bigger ball.
"""
function is_intersection_empty(B1::Ball2{N},
                               B2::Ball2{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool,Vector{N}}} where {N<:Real}
    center_diff_normed = norm(center(B2) - center(B1), 2)
    empty_intersection = center_diff_normed > B1.radius + B2.radius

    if !witness
        return empty_intersection
    elseif empty_intersection
        return (true, N[])
    end

    # compute a witness 'v' in the intersection
    if B1.radius <= B2.radius
        smaller = B1
        bigger = B2
    else
        smaller = B2
        bigger = B1
    end
    if center_diff_normed <= bigger.radius
        # smaller ball's center is contained in bigger ball
        v = smaller.center
    else
        # scale center difference with smaller ball's radius
	direction = (bigger.center - smaller.center)
        v = smaller.center + direction / center_diff_normed * smaller.radius
    end
    return (false, v)
end
