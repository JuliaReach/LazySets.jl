export Complement

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

isoperationtype(::Type{<:Complement}) = true

# special cases for which the complement is always convex
isconvextype(::Type{<:Complement{N,<:Union{EmptySet,HalfSpace,Universe}}}) where {N} = true

ispolyhedral(::Complement) = false
ispolyhedral(::Complement{N,<:Union{EmptySet,HalfSpace}}) where {N} = true

# the complement of the complement is the original set again
Complement(C::Complement) = C.X

"""
    dim(C::Complement)

Return the dimension of the complement of a set.

### Input

- `C` -- complement of a set

### Output

The ambient dimension of the complement of a set.
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
@validate function ∈(x::AbstractVector, C::Complement)
    return x ∉ C.X
end

"""
    isempty(C::Complement)

Check whether the complement of a set is empty.

### Input

- `C` -- complement of a set

### Output

`false` unless the original set is universal.

### Algorithm

We use the `isuniversal` function.
"""
function isempty(C::Complement)
    return isuniversal(C.X)
end

function isboundedtype(::Type{<:Complement{<:Real,<:Universe}})
    return true
end

function isboundedtype(::Type{<:Complement})
    return false
end

"""
    constraints_list(C::Complement)

Return the list of constraints of the complement of a set.

### Input

- `C` -- complement of a set

### Output

A vector of linear constraints.

### Notes

The method requires that the list of constraints of the complemented set can be
obtained. Then, each constraint is complemented and returned in the output
vector. The set union of this array corresponds to the concrete set complement.
"""
function constraints_list(C::Complement)
    @assert isconvextype(typeof(C.X)) "the constraints list of a complement " *
                                      "is only available for the complement of a convex polyhedron"

    clist = constraints_list(C.X)
    out = similar(clist)
    @inbounds for (i, ci) in enumerate(clist)
        out[i] = complement(ci)
    end
    return out
end

function translate(C::Complement, x::AbstractVector)
    return Complement(translate(C.X, x))
end

function concretize(C::Complement)
    return complement(concretize(C.X))
end
