const AbstractArraySet = Union{CartesianProductArray,
                               ConvexHullArray,
                               IntersectionArray,
                               MinkowskiSumArray,
                               UnionSetArray
                              }

function Base.getindex(X::AbstractArraySet, i)
    return getindex(array(X), i)
end

function Base.length(X::AbstractArraySet)
    return length(array(X))
end
