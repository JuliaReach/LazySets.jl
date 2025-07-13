"""
# Extended help

    linear_map(M::AbstractMatrix, X::Interval)

### Output

Either an interval or a zonotope, depending on the leading dimension (i.e., the
number of rows) of `M`:

- If `size(M, 1) == 1`, the output is an `Interval` obtained by scaling `X` by
  the matrix `M`.
- If `size(M, 1) ≠ 1`, the output is a `Zonotope` with center `M * center(X)`
  and the single generator `M * g`, where `g = (high(X)-low(X))/2`.
"""
@validate function linear_map(M::AbstractMatrix, X::Interval)
    nout = size(M, 1)
    if nout == 1
        return _linear_map_interval(M, X)
    else
        return _linear_map_zonotope(M, X)
    end
end

function _linear_map_interval(M::AbstractMatrix, X::Interval)
    α = @inbounds M[1, 1]
    return Interval(α * X.dat)
end

function _linear_map_zonotope(M::AbstractMatrix, X::Interval)
    require(@__MODULE__, :LazySets; fun_name="linear_map")

    nout = size(M, 1)
    cx = _center(X)
    gx = cx - min(X)
    N = promote_type(eltype(M), eltype(X))
    c = Vector{N}(undef, nout)
    gen = Matrix{N}(undef, nout, 1)
    @inbounds for i in 1:nout
        c[i] = M[i, 1] * cx
        gen[i] = M[i, 1] * gx
    end
    return Zonotope(c, gen)
end
