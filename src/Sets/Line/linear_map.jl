"""
# Extended help

    linear_map(M::AbstractMatrix, L::Line)

### Output

The line obtained by applying the linear map, if that still results in a line,
or a `Singleton` otherwise.

### Algorithm

We apply the linear map to the point and direction of `L`.
If the resulting direction is zero, the result is a singleton.
"""
@validate function linear_map(M::AbstractMatrix, L::Line)

    Mp = M * L.p
    Md = M * L.d
    if iszero(Md)
        require(@__MODULE__, :LazySets; fun_name="linear_map")

        return Singleton(Mp)
    end
    return Line(Mp, Md)
end
