import Base: isempty, ∈, ∪, Union

export UnionSet,
       UnionSetArray

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

# convenience constructor without type parameter
UnionSet(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} = UnionSet{N, S1, S2}(X, Y)

"""
    ∪

Alias for `UnionSet`.
"""
∪(X::LazySet, Y::LazySet) = UnionSet(X, Y)

"""
    dim(cup::Union)::Int

Return the dimension of the set union of two convex sets.

### Input

- `cup` -- union of two convex sets

### Output

The ambient dimension of the union of two convex sets.
"""
function dim(cup::Union)::Int
    return dim(cup.X)
end

"""
    UnionSetArray{N<:Real, S<:LazySet{N}}

Type that represents the set union of a finite number of convex sets.

### Fields

- `array` -- array of convex sets
"""
struct UnionSetArray{N<:Real, S<:LazySet{N}}
    array::Vector{S}

    # default constructor with dimension check
    function UnionSetArray{N, S}(X::Vector{S}) where {N<:Real, S<:LazySet{N}}
        @assert !isempty(X) "the array in this `UnionSetArray` is empty"
        n = dim(X[1])
        @assert all(x -> dim(x) == n, X) "sets in a union must have the same dimension"

        return new{N, S}(X)
    end
end

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
