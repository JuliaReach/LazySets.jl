import Base: ∈,
             isempty

export Complement,
       complement

"""
    Complement{N, S<:LazySet{N}}

Type that represents the complement of a set, that is the set

```math
Y = \\{y ∈ \\mathbb{R}^n : y ∉ X\\},
```
and it is often denoted with the ``C`` superscript, ``Y = X^C``.

### Fields

- `X` -- set

### Notes

Since `X` is assumed to be closed, unless `X` is empty or the universe, its
complement is open (i.e., not closed).
If `X` is empty, the universe, or a half-space, its complement is convex.

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
struct Complement{N, S<:LazySet{N}}
    X::S
end

isoperationtype(::Type{<:Complement}) = true

# the set complement is not convex in general
isconvextype(::Type{<:Complement}) = false

# special cases which are always convex
isconvextype(::Type{<:Complement{N, <:Union{EmptySet, HalfSpace, Universe}}}) where {N} = true

is_polyhedral(::Complement) where {N} = false
is_polyhedral(::Complement{N, <:Union{EmptySet, HalfSpace}}) where {N} = true

# the complement of the complement is the original set again
Complement(C::Complement) = C.X

"""
    dim(C::Complement)

Return the dimension of the complement of a set.

### Input

- `C` -- complement of a set

### Output

The dimension of the complement of a set.
"""
function dim(C::Complement)
    return dim(C.X)
end

"""
    ∈(x::AbstractVector, C::Complement)

Check whether a given point is contained in the complement of a set.

### Input

- `x` -- point/vector
- `C` -- complement of a set

### Output

`true` iff the vector is contained in the complement.

### Algorithm

```math
    x ∈ X^C ⟺ x ∉ X
```
"""
function ∈(x::AbstractVector, C::Complement)
    @assert length(x) == dim(C)
    return x ∉ C.X
end

"""
    isempty(C::Complement)

Return if the complement of a set is empty or not.

### Input

- `C` -- complement of a set

### Output

`false` unless the original set is universal.

### Algorithm

We use the `isuniversal` method.
"""
function isempty(C::Complement)
    return isuniversal(C.X)
end

function isboundedtype(::Type{<:Complement{<:Real, <:Universe}})
    return true
end

function isboundedtype(::Type{<:Complement})
    return false
end

# --  Fallback implementation, requires constraints list of C.X --

"""
    constraints_list(C::Complement)

Return the list of constraints of the complement of a set.

### Input

- `C` -- lazy set complement

### Output

A vector of linear constraints.

### Notes

The method requires that the list of constraints of the complemented set can be
obtained. Then, each constraint is complemented and returned in the output array.
The set union of such array corresponds to the concrete set complement.
"""
function constraints_list(C::Complement)
    clist = constraints_list(C.X)
    out = similar(clist)
    @inbounds for (i, ci) in enumerate(clist)
        out[i] = complement(ci)
    end
    return out
end
