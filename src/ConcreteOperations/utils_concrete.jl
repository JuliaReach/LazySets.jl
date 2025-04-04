"""
_concretize_lazy_array(arr, operation) -> LazySet
Helper to concretize lazy array operations (MinkowskiSumArray, CartesianProductArray).
"""
function _concretize_lazy_array(arr::AbstractVector, operation::Function)
    @assert !isempty(arr) "operation not supported on an empty array"
    x = first(arr)
    x = concretize(x)
    @inbounds for Y in @view arr[2:end]
        x = operation(x, concretize(Y))
    end
    return x
end