const AbstractArraySet = Union{CartesianProductArray,
                               ConvexHullArray,
                               IntersectionArray,
                               MinkowskiSumArray,
                               UnionSetArray
                              }

function Base.getindex(X::AbstractArraySet, i)
    return getindex(array(X), i)
end
