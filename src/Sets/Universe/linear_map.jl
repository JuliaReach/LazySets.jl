"""
# Extended help

    linear_map(M::AbstractMatrix, U::Universe)

### Output

Either a `Universe`, or an `HPolyhedron` if the map contains zero rows.
"""
@validate function linear_map(M::AbstractMatrix, U::Universe)
    m = size(M, 1)
    N = promote_type(eltype(M), eltype(U))
    first_zero_row = 0
    for i in 1:m
        if iszero(M[i, :])
            first_zero_row = i
            break
        end
    end
    if first_zero_row == 0
        return Universe{N}(m)
    end

    require(@__MODULE__, :LazySets; fun_name="linear_map")
    clist = Vector{HalfSpace{N,SingleEntryVector{N}}}()
    _push_zero_constraint!(clist, N, first_zero_row, m)
    for i in (first_zero_row + 1):m
        if iszero(M[i, :])
            _push_zero_constraint!(clist, N, i, m)
        end
    end
    return HPolyhedron(clist)
end

function _push_zero_constraint!(clist, N, i, n)
    push!(clist, HalfSpace(SingleEntryVector(i, n, one(N)), zero(N)))
    push!(clist, HalfSpace(SingleEntryVector(i, n, N(-1)), zero(N)))
end
