"""
    Line{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a line of the form

```math
    \\{y ∈ ℝ^n: y = p + λd, λ ∈ ℝ\\}
```
where ``p`` is a point on the line and ``d`` is its direction vector (not
necessarily normalized).

### Fields

- `p` -- point on the line
- `d` -- direction

### Examples

There are three constructors. The optional keyword argument `normalize`
(default: `false`) can be used to normalize the direction of the resulting line
to have norm 1 (w.r.t. the Euclidean norm).

1) The default constructor takes the fields `p` and `d`:

The line passing through the point ``[-1, 2, 3]`` and parallel to the vector
``[3, 0, -1]``:

```jldoctest
julia> Line([-1.0, 2, 3], [3.0, 0, -1])
Line{Float64, Vector{Float64}}([-1.0, 2.0, 3.0], [3.0, 0.0, -1.0])

julia> Line([-1.0, 2, 3], [3.0, 0, -1]; normalize=true)
Line{Float64, Vector{Float64}}([-1.0, 2.0, 3.0], [0.9486832980505138, 0.0, -0.31622776601683794])
```

2) The second constructor takes two points, `from` and `to`, as keyword
arguments, and returns the line through them. See the algorithm section for
details.

```jldoctest
julia> Line(from=[-1.0, 2, 3], to=[-4.0, 2, 4])
Line{Float64, Vector{Float64}}([-1.0, 2.0, 3.0], [3.0, 0.0, -1.0])
```

3) The third constructor resembles `Line2D` and only works for two-dimensional
lines. It takes two inputs, `a` and `b`, and constructs the line such that
``a ⋅ x = b``.

```jldoctest
julia> Line([2.0, 0], 1.)
Line{Float64, Vector{Float64}}([0.5, 0.0], [0.0, -1.0])
```

### Algorithm

Given two points ``p ∈ ℝ^n`` and ``q ∈ ℝ^n``, the line that
passes through these two points is
`L: `\\{y ∈ ℝ^n: y = p + λ(q - p), λ ∈ ℝ\\}``.
"""
struct Line{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    p::VN
    d::VN

    # default constructor with length constraint
    function Line(p::VN, d::VN; check_direction::Bool=true,
                  normalize=false) where {N,VN<:AbstractVector{N}}
        if check_direction
            @assert !iszero(d) "a line needs a non-zero direction vector"
        end
        d_n = normalize ? LinearAlgebra.normalize(d) : d
        return new{N,VN}(p, d_n)
    end
end

function Line(; from::AbstractVector, to::AbstractVector, normalize=false)
    d = from - to
    @assert !iszero(d) "points `$from` and `$to` should be distinct"
    return Line(from, d; normalize=normalize)
end

function Line(a::AbstractVector{N}, b::N; normalize=false) where {N}
    @assert length(a) == 2 "expected a normal vector of length two, but it " *
                           "is $(length(a))-dimensional"

    got_horizontal = iszero(a[1])
    got_vertical = iszero(a[2])

    if got_horizontal && got_vertical
        throw(ArgumentError("the vector $a must be non-zero"))
    end

    if got_horizontal
        α = b / a[2]
        p = [zero(N), α]
        q = [one(N), α]
    elseif got_vertical
        β = b / a[1]
        p = [β, zero(N)]
        q = [β, one(N)]
    else
        α = b / a[2]
        μ = a[1] / a[2]
        p = [zero(N), α]
        q = [one(N), α - μ]
    end
    return Line(; from=p, to=q, normalize=normalize)
end
