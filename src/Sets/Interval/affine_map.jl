function affine_map(M, X::Interval, v::AbstractVector; kwargs...)
    @assert size(M, 2) == 1 "cannot apply an affine map of dimension $(size(M, 2)) " *
                            "to an interval"
    @assert size(M, 1) == length(v) "cannot apply an affine map of matrix dimension " *
                                    "$(size(M, 1)) and translation dimension $(length(v))"
    @inbounds if length(v) == 1
        return Interval(M[1, 1] * X.dat + v[1])
    end
    return translate(linear_map(M, X; kwargs...), v)
end
