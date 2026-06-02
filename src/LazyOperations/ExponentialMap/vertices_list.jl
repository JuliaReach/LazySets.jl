"""
    vertices_list(em::ExponentialMap; [backend]=get_exponential_backend())

Return the list of vertices of a (polytopic) exponential map.

### Input

- `em`      -- polytopic exponential map
- `backend` -- (optional; default: `get_exponential_backend()`) exponentiation
               backend

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polytopic.
Then the result is just the exponential map applied to the vertices of `X`.
"""
@validate function vertices_list(em::ExponentialMap; backend=get_exponential_backend())
    # collect vertices lists of wrapped set
    vlist_X = vertices_list(em.X)

    # create resulting vertices list
    N = eltype(em)
    vlist = Vector{Vector{N}}(undef, length(vlist_X))
    @inbounds for (i, v) in enumerate(vlist_X)
        vlist[i] = _expmv(backend, one(N), em.expmat.M, v)
    end

    return vlist
end
