"""
# Extended help

    linear_map(M::AbstractMatrix, x::Interval)

### Output

Either an interval or a zonotope, depending on the leading dimension (i.e., the
number of rows) of `M`:

- If `size(M, 1) == 1`, the output is an `Interval` obtained by scaling `x` by
  the matrix `M`.
- If `size(M, 1) ≠ 1`, the output is a `Zonotope` with center `M * center(x)`
  and the single generator `M * g`, where `g = (high(x)-low(x))/2`.
"""
function linear_map(M::AbstractMatrix, x::Interval)
    @assert size(M, 2) == 1 "a linear map of size $(size(M)) " *
                            "cannot be applied to an interval"
    nout = size(M, 1)
    if nout == 1
        return _linear_map_interval(M, x)
    else
        return _linear_map_zonotope(M, x)
    end
end

function _linear_map_interval(M::AbstractMatrix, x::Interval)
    α = @inbounds M[1, 1]
    return Interval(α * x.dat)
end

function _linear_map_zonotope(M::AbstractMatrix, x::Interval)
    require(@__MODULE__, :LazySets; fun_name="linear_map")

    nout = size(M, 1)
    cx = _center(x)
    gx = cx - min(x)
    N = promote_type(eltype(M), eltype(x))
    c = Vector{N}(undef, nout)
    gen = Matrix{N}(undef, nout, 1)
    @inbounds for i in 1:nout
        c[i] = M[i, 1] * cx
        gen[i] = M[i, 1] * gx
    end
    return Zonotope(c, gen)
end
