"""
# Extended help

    in(x::AbstractVector, H::Hyperplane)

### Algorithm

We just check whether ``x`` satisfies ``a⋅x = b``.
"""
@validate function in(x::AbstractVector, H::Hyperplane)
    return _isapprox(dot(H.a, x), H.b)
end
