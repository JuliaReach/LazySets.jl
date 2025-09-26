@validate function constraints_list(P::HParallelotope)
    D, c = P.directions, P.offset
    N, VN = _parameters(P)
    return _constraints_list_hparallelotope(D, c, N, VN)
end

# reason: `Documenter` cannot annotate `constraints_list` with type parameters
function _parameters(::HParallelotope{N,VN}) where {N,VN}
    return (N, VN)
end

function _constraints_list_hparallelotope(D, c, N, VN)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    if isempty(D)
        return Vector{HalfSpace{N,VN}}(undef, 0)
    end
    n = size(D, 1)
    clist = Vector{HalfSpace{N,VN}}(undef, 2n)
    @inbounds for i in 1:n
        clist[i] = HalfSpace(D[i, :], c[i])
        clist[i + n] = HalfSpace(-D[i, :], c[i + n])
    end
    return clist
end
