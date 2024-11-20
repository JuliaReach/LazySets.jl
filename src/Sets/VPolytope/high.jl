function high(P::VPolytope)
    return _high_vlist(P)
end

function high(P::VPolytope, i::Int)
    @assert 1 <= i <= dim(P) "invalid index $i for set of dimension $(dim(P))"
    return _high_vlist(P, i)
end
