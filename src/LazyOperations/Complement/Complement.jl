"""
    Complement{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the complement of a set, i.e., the set

```math
Y = \\{y ∈ ℝ^n : y ∉ X\\}.
```

The complement is often denoted with the ``C`` superscript, as in ``Y = X^C``.

### Fields

- `X` -- set

### Notes

If `X` is empty, the universe, or a half-space, its complement is convex.

Since `X` is assumed to be closed, unless `X` is empty or the universe, its
complement is open (i.e., not closed). In this library, all sets are closed, so
the set is usually not represented exactly at the boundary.

The complement of the complement is the original set again.

### Examples

```jldoctest
julia> B = BallInf(zeros(2), 1.);

julia> C = Complement(B)
Complement{Float64, BallInf{Float64, Vector{Float64}}}(BallInf{Float64, Vector{Float64}}([0.0, 0.0], 1.0))

julia> Complement(C)
BallInf{Float64, Vector{Float64}}([0.0, 0.0], 1.0)
```
"""
struct Complement{N,S<:LazySet{N}} <: LazySet{N}
    X::S
end

# the complement of the complement is the original set again
Complement(C::Complement) = C.X

include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedraltype.jl")
include("in.jl")
include("translate.jl")
