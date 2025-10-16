import Base: +, getindex

export MinkowskiSum, ⊕,
       MinkowskiSum!,
       swap

"""
    MinkowskiSum{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of two sets, i.e., the set

```math
X ⊕ Y = \\{x + y : x ∈ X, y ∈ Y\\}.
```

### Fields

- `X` -- set
- `Y` -- set

### Notes

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSum`.

The Minkowski sum preserves convexity: if the set arguments are convex, then
their Minkowski sum is convex as well.

The convenience aliases `⊕` and `+` are also available. `⊕` can be typed by
`\\oplus<tab>`.
"""
struct MinkowskiSum{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function MinkowskiSum(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a Minkowski sum must have the same dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

+(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

⊕(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

isoperationtype(::Type{<:MinkowskiSum}) = true
concrete_function(::Type{<:MinkowskiSum}) = minkowski_sum

isconvextype(::Type{MinkowskiSum{N,S1,S2}}) where {N,S1,S2} = isconvextype(S1) && isconvextype(S2)

function ispolyhedraltype(::Type{<:MinkowskiSum{N,S1,S2}}) where {N,S1,S2}
    return ispolyhedraltype(S1) && ispolyhedraltype(S2)
end

function ispolyhedral(ms::MinkowskiSum)
    ispolyhedraltype(typeof(ms)) && return true

    return ispolyhedral(ms.X) && ispolyhedral(ms.Y)
end

# ZeroSet is the neutral element for MinkowskiSum
@neutral(MinkowskiSum, ZeroSet)

# EmptySet and Universe are the absorbing elements for MinkowskiSum
@absorbing(MinkowskiSum, EmptySet)
# @absorbing(MinkowskiSum, Universe)  # TODO problematic

# interface for binary set operations
Base.first(ms::MinkowskiSum) = ms.X
second(ms::MinkowskiSum) = ms.Y
@declare_binary_operation(MinkowskiSum)

"""
    swap(ms::MinkowskiSum)

Return a new `MinkowskiSum` object with the arguments swapped.

### Input

- `ms` -- Minkowski sum of two sets

### Output

A new `MinkowskiSum` object with the arguments swapped.
"""
function swap(ms::MinkowskiSum)
    return MinkowskiSum(ms.Y, ms.X)
end

"""
    dim(ms::MinkowskiSum)

Return the dimension of a Minkowski sum of two sets.

### Input

- `ms` -- Minkowski sum of two sets

### Output

The ambient dimension of the Minkowski sum of two sets.
"""
function dim(ms::MinkowskiSum)
    return dim(ms.X)
end

"""
    σ(d::AbstractVector, ms::MinkowskiSum)

Return a support vector of a Minkowski sum of two sets.

### Input

- `d`  -- direction
- `ms` -- Minkowski sum of two sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.

### Algorithm

A valid support vector in direction ``d`` of the Minkowski sum of two sets ``X``
and ``Y`` is the sum of the support vectors of ``X`` and ``Y`` in direction
``d``.
"""
@validate function σ(d::AbstractVector, ms::MinkowskiSum)
    return σ(d, ms.X) + σ(d, ms.Y)
end

"""
    ρ(d::AbstractVector, ms::MinkowskiSum)

Evaluate the support function of a Minkowski sum of two sets.

### Input

- `d`  -- direction
- `ms` -- Minkowski sum of two sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function in direction ``d`` of the Minkowski sum of two sets ``X``
and ``Y`` is the sum of the support functions of ``X`` and ``Y`` in direction
``d``.
"""
@validate function ρ(d::AbstractVector, ms::MinkowskiSum)
    return ρ(d, ms.X) + ρ(d, ms.Y)
end

"""
    isbounded(ms::MinkowskiSum)

Check whether a Minkowski sum of two sets is bounded.

### Input

- `ms` -- Minkowski sum of two sets

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(ms::MinkowskiSum)
    return isbounded(ms.X) && isbounded(ms.Y)
end

function isboundedtype(::Type{MinkowskiSum{N,S1,S2}}) where {N,S1,S2}
    return isboundedtype(S1) && isboundedtype(S2)
end

"""
    isempty(ms::MinkowskiSum)

Check whether a Minkowski sum of two sets is empty.

### Input

- `ms` -- Minkowski sum of two sets

### Output

`true` iff any of the wrapped sets are empty.
"""
function isempty(ms::MinkowskiSum)
    return isempty(ms.X) || isempty(ms.Y)
end

"""
    center(ms::MinkowskiSum)

Return the center of a Minkowski sum of two centrally-symmetric sets.

### Input

- `ms` -- Minkowski sum of two centrally-symmetric sets

### Output

The center of the Minkowski sum.
"""
function center(ms::MinkowskiSum)
    return center(ms.X) + center(ms.Y)
end

"""
    constraints_list(ms::MinkowskiSum)

Return a list of constraints of the Minkowski sum of two polyhedral sets.

### Input

- `ms` -- Minkowski sum of two polyhedral sets

### Output

The list of constraints of the Minkowski sum.

### Algorithm

We compute a concrete set representation via `minkowski_sum` and call
`constraints_list` on the result.
"""
function constraints_list(ms::MinkowskiSum)
    return constraints_list(minkowski_sum(ms.X, ms.Y))
end

"""
    ∈(x::AbstractVector, ms::MinkowskiSum{N,<:AbstractSingleton}) where {N}

Check whether a given point is contained in the Minkowski sum of a singleton
and another set.

### Input

- `x`  -- point/vector
- `ms` -- Minkowski sum of a singleton and another set

### Output

`true` iff ``x ∈ ms``.

### Algorithm

Note that ``x ∈ (S ⊕ P)``, where ``S = \\{s\\}``  is a singleton set and
``P`` is a set, if and only if ``(x-s) ∈ P``.
"""
@validate function ∈(x::AbstractVector, ms::MinkowskiSum{N,<:AbstractSingleton}) where {N}
    return _in_singleton_msum(x, ms.X, ms.Y)
end

# symmetric method
@validate function ∈(x::AbstractVector,
                     ms::MinkowskiSum{N,<:LazySet,<:AbstractSingleton}) where {N}
    return _in_singleton_msum(x, ms.Y, ms.X)
end

# disambiguation
@validate function ∈(x::AbstractVector,
                     ms::MinkowskiSum{N,<:AbstractSingleton,<:AbstractSingleton}) where {N}
    return _in_singleton_msum(x, ms.X, ms.Y)
end

@inline _in_singleton_msum(x, X, Y) = (x - element(X)) ∈ Y

"""
    vertices_list(ms::MinkowskiSum)

Return a list of vertices for the Minkowski sum of two sets.

### Input

- `ms` -- Minkowski sum of two sets

### Output

A list of vertices of the Minkowski sum of two sets.

### Algorithm

We compute the concrete Minkowski sum (via `minkowski_sum`) and call
`vertices_list` on the result.
"""
function vertices_list(ms::MinkowskiSum)
    return vertices_list(minkowski_sum(ms.X, ms.Y))
end

@validate function translate(ms::MinkowskiSum, x::AbstractVector)
    X = translate(first(ms), x)
    Y = translate(second(ms), x)
    return MinkowskiSum(X, Y)
end
