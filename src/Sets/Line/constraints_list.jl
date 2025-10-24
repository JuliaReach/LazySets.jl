"""
# Extended help

    constraints_list(L::Line)

### Output

A list containing `2n-2` half-spaces whose intersection is `L`, where `n` is the
ambient dimension of `L`.
"""
function constraints_list(L::Line)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    p = L.p
    n = length(p)
    d = reshape(L.d, 1, n)
    K = nullspace(d)
    m = size(K, 2)
    @assert m == n - 1 "expected $(n - 1) normal half-spaces, got $m"

    N, VN = _parameters(L)
    out = Vector{HalfSpace{N,VN}}(undef, 2m)
    idx = 1
    @inbounds for Kj in eachcol(K)
        b = dot(Kj, p)
        out[idx] = HalfSpace(Kj, b)
        out[idx + 1] = HalfSpace(-Kj, -b)
        idx += 2
    end
    return out
end

# reason: Documenter cannot annotate `constraints_list` with type parameters
function _parameters(::Line{N,VN}) where {N,VN}
    return (N, VN)
end
