import Base: isempty, ∈, ∪

export UnionSet,
       UnionSetArray,
       array,
       swap

# ========================================
# Binary set union
# ========================================

"""
    UnionSet{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}

Type that represents the set union of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set
"""
struct UnionSet{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
    X::S1
    Y::S2

    # default constructor with dimension check
    function UnionSet{N, S1, S2}(X::S1, Y::S2) where{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in a union must have the same dimension"
        return new{N, S1, S2}(X, Y)
    end
end

isoperationtype(::Type{<:UnionSet}) = true

# convenience constructor without type parameter
UnionSet(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} = UnionSet{N, S1, S2}(X, Y)

# EmptySet is the neutral element for UnionSet
@neutral(UnionSet, EmptySet)

# Universe is the absorbing element for UnionSet
@absorbing(UnionSet, Universe)

"""
    ∪

Alias for `UnionSet`.
"""
∪(X::LazySet, Y::LazySet) = UnionSet(X, Y)

"""
    swap(cup::UnionSet)

Return a new `UnionSet` object with the arguments swapped.

### Input

- `cup` -- union of two convex sets

### Output

A new `UnionSet` object with the arguments swapped.
"""
function swap(cup::UnionSet)
    return UnionSet(cup.Y, cup.X)
end

"""
    dim(cup::UnionSet)::Int

Return the dimension of the set union of two convex sets.

### Input

- `cup` -- union of two convex sets

### Output

The ambient dimension of the union of two convex sets.
"""
function dim(cup::UnionSet)::Int
    return dim(cup.X)
end

"""
    σ(d::AbstractVector{N}, cup::UnionSet{N}; [algorithm]="support_vector") where {N<:Real}

Return the support vector of the union of two convex sets in a given direction.

### Input

- `d`         -- direction
- `cup`       -- union of two convex sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                 the support vector; if "support_vector", use the support
                 vector of each argument; if "support_function" use the support
                 function of each argument and evaluate the support vector of only
                 one of them

### Output

The support vector in the given direction.

### Algorithm

The support vector of the union of two convex sets ``X`` and ``Y`` can be obtained
as the vector that maximizes the support function of either ``X`` or ``Y``, i.e.
it is sufficient to find the ``\\argmax(ρ(d, X), ρ(d, Y)])`` and evaluate its support
vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of ``X`` and ``Y`` and then compares the support function using
a dot product. If it happens that the support function can be more efficiently
computed (without passing through the support vector), consider using the alternative
`algorithm="support_function"` implementation, which evaluates the support function
of each set directly and then calls only the support vector of either ``X`` *or*
``Y``.
"""
function σ(d::AbstractVector{N}, cup::UnionSet{N}; algorithm="support_vector") where {N<:Real}
    X, Y = cup.X, cup.Y
    if algorithm == "support_vector"
        σX, σY = σ(d, X), σ(d, Y)
        return dot(d, σX) > dot(d, σY) ? σX : σY
    elseif algorithm == "support_function"
        m = argmax([ρ(d, X), ρ(d, Y)])
        return m == 1 ? σ(d, X) : σ(d, Y)
    else
        error("algorithm $algorithm for the support vector of a `UnionSet` is unknown")
    end
end

"""
    ρ(d::AbstractVector{N}, cup::UnionSet{N}) where {N<:Real}

Return the support function of the union of two convex sets in a given direction.

### Input

- `d`   -- direction
- `cup` -- union of two convex sets

### Output

The support function in the given direction.

### Algorithm

The support function of the union of two convex sets ``X`` and ``Y`` is the
maximum of the support functions of ``X`` and ``Y``.
"""
function ρ(d::AbstractVector{N}, cup::UnionSet{N}) where {N<:Real}
    X, Y = cup.X, cup.Y
    return max(ρ(d, X), ρ(d, Y))
end

"""
    an_element(cup::UnionSet{N})::Vector{N} where {N<:Real}

Return some element of a union of two convex sets.

### Input

- `cup` -- union of two convex sets

### Output

An element in the union of two convex sets.

### Algorithm

We use `an_element` on the first wrapped set.
"""
function an_element(cup::UnionSet{N})::Vector{N} where {N<:Real}
    return an_element(cup.X)
end

"""
    ∈(x::AbstractVector{N}, cup::UnionSet{N})::Bool where {N<:Real}

Check whether a given point is contained in a union of two convex sets.

### Input

- `x`   -- point/vector
- `cup` -- union of two convex sets

### Output

`true` iff ``x ∈ cup``.
"""
function ∈(x::AbstractVector{N}, cup::UnionSet{N})::Bool where {N<:Real}
    return x ∈ cup.X || x ∈ cup.Y
end

"""
    isempty(cup::UnionSet)::Bool

Check whether a union of two convex sets is empty.

### Input

- `cup` -- union of two convex sets

### Output

`true` iff the union is empty.
"""
function isempty(cup::UnionSet)::Bool
    return isempty(cup.X) && isempty(cup.Y)
end

"""
    isbounded(cup::UnionSet)::Bool

Determine whether a union of two convex sets is bounded.

### Input

- `cup` -- union of two convex sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSet)::Bool
    return isbounded(cup.X) && isbounded(cup.Y)
end

# ========================================
# n-ary set union
# ========================================

"""
    UnionSetArray{N<:Real, S<:LazySet{N}}

Type that represents the set union of a finite number of convex sets.

### Fields

- `array` -- array of convex sets
"""
struct UnionSetArray{N<:Real, S<:LazySet{N}}
    array::Vector{S}
end

isoperationtype(::Type{<:UnionSetArray}) = true

# EmptySet is the neutral element for UnionSetArray
@neutral(UnionSetArray, EmptySet)

# Universe is the absorbing element for UnionSetArray
@absorbing(UnionSetArray, Universe)

# add functions connecting UnionSet and UnionSetArray
@declare_array_version(UnionSet, UnionSetArray)

"""
    dim(cup::UnionSetArray)::Int

Return the dimension of the set union of a finite number of convex sets.

### Input

- `cup` -- union of a finite number of convex sets

### Output

The ambient dimension of the union of a finite number of convex sets.
"""
function dim(cup::UnionSetArray)::Int
    return dim(cup.array[1])
end

"""
    array(cup::UnionSetArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a union of a finite number of convex sets.

### Input

- `cup` -- union of a finite number of convex sets

### Output

The array that holds the union of a finite number of convex sets.
"""
function array(cup::UnionSetArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return cup.array
end

"""
    σ(d::AbstractVector{N}, cup::UnionSetArray{N}; [algorithm]="support_vector") where {N<:Real}

Return the support vector of the union of a finite number of convex sets in
a given direction.

### Input

- `d`         -- direction
- `cup`       -- union of a finite number of convex sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                 the support vector; if "support_vector", use the support
                 vector of each argument; if "support_function" use the support
                 function of each argument and evaluate the support vector of only
                 one of them

### Output

The support vector in the given direction.

### Algorithm

The support vector of the union of a finite number of convex sets ``X₁, X₂, ...``
can be obtained as the vector that maximizes the support function, i.e.
it is sufficient to find the ``\\argmax(ρ(d, X₂), ρ(d, X₂), ...])`` and evaluate
its support vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of all ``X₁, X₂, ...`` and then compares the support function using
a dot product. If it happens that the support function can be more efficiently
computed (without passing through the support vector), consider using the alternative
`algorithm="support_function"` implementation, which evaluates the support function
of each set directly and then calls only the support vector of one of the ``Xᵢ``.
"""
function σ(d::AbstractVector{N}, cup::UnionSetArray{N}; algorithm="support_vector") where {N<:Real}
    A = array(cup)
    if algorithm == "support_vector"
        σarray = map(Xi -> σ(d, Xi), A)
        ρarray = map(vi -> dot(d, vi), σarray)
        m = argmax(ρarray)
        return σarray[m]
    elseif algorithm == "support_function"
        ρarray = map(Xi -> ρ(d, Xi), A)
        m = argmax(ρarray)
        return σ(d, A[m])
    else
        error("algorithm $algorithm for the support vector of a `UnionSetArray` is unknown")
    end
end

"""
    ρ(d::AbstractVector{N}, cup::UnionSetArray{N}) where {N<:Real}

Return the support function of the union of a finite number of convex sets in
a given direction.

### Input

- `d`   -- direction
- `cup` -- union of a finite number of convex sets

### Output

The support function in the given direction.

### Algorithm

The support function of the union of a finite number of convex sets ``X₁, X₂, ...``
can be obtained as the maximum of ``ρ(d, X₂), ρ(d, X₂), ...``.
"""
function ρ(d::AbstractVector{N}, cup::UnionSetArray{N}) where {N<:Real}
    A = array(cup)
    ρarray = map(Xi -> ρ(d, Xi), A)
    return maximum(ρarray)
end

"""
    an_element(cup::UnionSetArray{N})::Vector{N} where {N<:Real}

Return some element of a union of a finite number of convex sets.

### Input

- `cup` -- union of a finite number of convex sets

### Output

An element in the union of a finite number of convex sets.

### Algorithm

We use `an_element` on the first wrapped set.
"""
function an_element(cup::UnionSetArray{N})::Vector{N} where {N<:Real}
    return an_element(array(cup)[1])
end

"""
    ∈(x::AbstractVector{N}, cup::UnionSetArray{N})::Bool where {N<:Real}

Check whether a given point is contained in a union of a finite number of convex
sets.

### Input

- `x`   -- point/vector
- `cup` -- union of a finite number of convex sets

### Output

`true` iff ``x ∈ cup``.
"""
function ∈(x::AbstractVector{N}, cup::UnionSetArray{N})::Bool where {N<:Real}
    return any(X -> x ∈ X, array(cup))
end

"""
    isempty(cup::UnionSetArray)::Bool

Check whether a union of a finite number of convex sets is empty.

### Input

- `cup` -- union of a finite number of convex sets

### Output

`true` iff the union is empty.
"""
function isempty(cup::UnionSetArray)::Bool
    return all(isempty, array(cup))
end

"""
    isbounded(cup::UnionSetArray)::Bool

Determine whether a union of a finite number of convex sets is bounded.

### Input

- `cup` -- union of a finite number of convex sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSetArray)::Bool
    return all(isbounded, array(cup))
end
