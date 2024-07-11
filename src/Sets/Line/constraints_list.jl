"""
    constraints_list(L::Line)

Return the list of constraints of a line.

### Input

- `L` -- line

### Output

A list containing `2n-2` half-spaces whose intersection is `L`, where `n` is the
ambient dimension of `L`.
"""
function constraints_list(L::Line)
    p = L.p
    n = length(p)
    d = reshape(L.d, 1, n)
    K = nullspace(d)
    m = size(K, 2)
    @assert m == n - 1 "expected $(n - 1) normal half-spaces, got $m"

    N, VN = _parameters(L)
    out = Vector{HalfSpace{N,VN}}(undef, 2m)
    idx = 1
    @inbounds for j in 1:m
        Kj = K[:, j]
        b = dot(Kj, p)
        out[idx] = HalfSpace(Kj, b)
        out[idx + 1] = HalfSpace(-Kj, -b)
        idx += 2
    end
    return out
end

function _parameters(::Line{N,VN}) where {N,VN}
    return (N, VN)
end
