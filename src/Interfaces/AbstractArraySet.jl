const AbstractArraySet = Union{CartesianProductArray,
                               ConvexHullArray,
                               IntersectionArray,
                               MinkowskiSumArray,
                               UnionSetArray}

function Base.getindex(X::AbstractArraySet, i::Int)
    return getindex(array(X), i)
end

function Base.getindex(X::AbstractArraySet, indices::AbstractVector{Int})
    return [X[i] for i in indices]
end

function Base.lastindex(X::AbstractArraySet)
    return length(array(X))
end

function Base.length(X::AbstractArraySet)
    return length(array(X))
end

function Base.iterate(X::AbstractArraySet, state=1)
    return iterate(array(X), state)
end

"""
    flatten(X::LazySets.AbstractArraySet)

Flatten an array set, i.e., resolve potential nestings.

### Examples

```jldoctest
julia> E1 = EmptySet(1); E2 = EmptySet(2); E3 = EmptySet(3);

julia> X = MinkowskiSumArray([E1, MinkowskiSumArray([E2, E2])])
MinkowskiSumArray{Float64, LazySet{Float64}}(LazySet{Float64}[∅(1), MinkowskiSumArray{Float64, ∅}(∅[∅(2), ∅(2)])])

julia> flatten(X)
MinkowskiSumArray{Float64, ∅}(∅[∅(1), ∅(2), ∅(2)])
```
"""
function flatten end

function flatten!(arr, X, bin_op)
    if X isa bin_op
        flatten!(arr, first(X), bin_op)
        flatten!(arr, second(X), bin_op)
    elseif X isa array_constructor(bin_op)
        for Xi in X
            flatten!(arr, Xi, bin_op)
        end
    else
        push!(arr, X)
    end
    return arr
end

function _concretize_lazy_array(container::T) where {T<:AbstractArraySet}
    arr = array(container)
    @assert !isempty(arr) "`concretize` not supported on empty array"

    op = concrete_function(T)
    X = concretize(first(arr))
    @inbounds for Y in @view arr[2:end]
        X = op(X, concretize(Y))
    end
    return X
end
