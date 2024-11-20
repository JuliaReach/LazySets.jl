"""
# Extended help

    minkowski_difference(I1::Interval, I2::Interval)

### Output

An `Interval`, or an `EmptySet` if the difference is empty.
"""
function minkowski_difference(I1::Interval, I2::Interval)
    l = min(I1) - min(I2)
    h = max(I1) - max(I2)
    if h < l
        require(@__MODULE__, :LazySets; fun_name="minkowski_difference")

        N = promote_type(eltype(I1), eltype(I2))
        return EmptySet{N}(1)
    end
    return Interval(l, h)
end
