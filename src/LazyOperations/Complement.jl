import Base: ∈,
             isempty

export Complement,
       complement

"""
    Complement{N, S<:LazySet{N}}

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
julia> B = BallInf(zeros(2), 1.);

julia> C = Complement(B)
Complement{Float64,BallInf{Float64,Array{Float64,1}}}(BallInf{Float64,Array{Float64,1}}([0.0, 0.0], 1.0))

julia> Complement(C)
BallInf{Float64,Array{Float64,1}}([0.0, 0.0], 1.0)
```
"""
struct Complement{N, S<:LazySet{N}}
    X::S
end

isoperationtype(::Type{<:Complement}) = true

# the set complement is not convex in general
isconvextype(::Type{<:Complement}) = false

# special cases which are always convex
isconvextype(::Type{Complement{N, HalfSpace{N, VN}}}) where {N, VN} = true
isconvextype(::Type{Complement{N, EmptySet{N}}}) where {N} = true
isconvextype(::Type{Complement{N, Universe{N}}}) where {N} = true

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
function dim(C::Complement)
    return dim(C.X)
end

"""
    ∈(x::AbstractVector, C::Complement)

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
function ∈(x::AbstractVector, C::Complement)
    @assert length(x) == dim(C)
    return x ∉ C.X
end

"""
    isempty(C::Complement)

Return if the complement of a convex set is empty or not.

### Input

- `C` -- complement of a convex set

### Output

`false` unless the original set is universal.

### Algorithm

We use the `isuniversal` method.
"""
function isempty(C::Complement)
    return isuniversal(C.X)
end

# -- Concrete complement implementations --

"""
    complement(X::LazySet)

Return the complement of a set.

### Input

- `X` -- set

### Output

A `UnionSetArray` of half-spaces, i.e. the output is the the union of the linear
constraints which are obtained by complementing each constraint of `X`.

### Algorithm

The principle used in this function is that if ``X`` and ``Y`` are any pair of sets,
then ``(X ∩ Y)^C = X^C ∪ Y^C``. In particular, we can apply this rule for each constraint
that defines a polyhedral set, hence the concrete complement can be represented as the set
union of the complement of each constraint.
"""
function complement(X::LazySet)
    return UnionSetArray(constraints_list(Complement(X)))
end

"""
    complement(H::HalfSpace)

Return the complement of a half-space.

### Input

- `H` -- halfspace

### Output

The halfspace that is complementary to `H`. If ``H: \\langle a, x \\rangle ≤ b``,
then this function returns the halfspace ``H′: \\langle a, x \\rangle ≥ b``.
(Note that complementarity is understood in a relaxed sense, since the intersection
of ``H`` and ``H′`` is non-empty).
"""
function complement(H::HalfSpace)
    return HalfSpace(-H.a, -H.b)
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
    for (i, ci) in enumerate(clist)
        out[i] = complement(ci)
    end
    return out
end
