import Base: ∈,
             isempty

export Complement

"""
    Complement{N<:Real, S<:LazySet{N}}

Type that represents the complement of a convex set.

### Fields

- `X` -- convex set

### Notes

Since `X` is assumed to be closed, unless `X` is empty or the universe, its
complement is open (i.e., not closed).
Since `X` is assumed to be closed, unless `X` is empty, the universe, or a
half-space, its complement is not convex.

The complement of the complement is the original set again.

### Examples

```jldoctest
julia> B = BallInf(zeros(2), 1.); C = Complemeisconvextype(::Type{<:Complement) = falsent(B)
Complement{Float64,BallInf{Float64}}(BallInf{Float64}([0.0, 0.0], 1.0))

julia> Complement(C)
BallInf{Float64}([0.0, 0.0], 1.0)
```
"""
struct Complement{N<:Real, S<:LazySet{N}}
    X::S
end

isoperationtype(::Type{<:Complement}) = true
isconvextype(::Type{<:Complement}) = false

# the complement of the complement is the original set again
Complement(C::Complement) = C.X

"""
    dim(C::Complement)

Return the dimension of the complement of a convex set.

### Input

- `C` -- complement of a convex set

### Output

The dimension of the complement of a convex set.
"""
function dim(C::Complement)::Int
    return dim(C.X)
end

"""
    ∈(x::AbstractVector{N}, C::Complement{N})::Bool where {N<:Real}

Check whether a given point is contained in the complement of a convex set.

### Input

- `x` -- point/vector
- `C` -- complement of a convex set

### Output

`true` iff the vector is contained in the complement.

### Algorithm

```math
    x ∈ X^C ⟺ x ∉ X
```
"""
function ∈(x::AbstractVector{N}, C::Complement{N})::Bool where {N<:Real}
    @assert length(x) == dim(C)
    return x ∉ C.X
end

"""
    isempty(C::Complement)::Bool

Return if the complement of a convex set is empty or not.

### Input

- `C` -- complement of a convex set

### Output

`false` unless the original set is universal.

### Algorithm

We use the `isuniversal` method.
"""
function isempty(C::Complement)::Bool
    return isuniversal(C.X)
end
