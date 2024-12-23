"""
# Extended help

    isdisjoint(H1::HalfSpace, H2::HalfSpace, [witness]::Bool=false)

### Algorithm

Two half-spaces do not intersect if and only if their normal vectors point in
the opposite direction and there is a gap between the two defining hyperplanes.

The latter can be checked as follows:
Let ``H1 : a_1⋅x = b_1`` and ``H2 : a_2⋅x = b_2``.
Then we already know that ``a_2 = -k⋅a_1`` for some positive scaling factor
``k``.
Let ``x_1`` be a point on the defining hyperplane of ``H1``.
We construct a line segment from ``x_1`` to the point ``x_2`` on the defining
hyperplane of ``hs_2`` by shooting a ray from ``x_1`` with direction ``a_1``.
Thus we look for a factor ``s`` such that ``(x_1 + s⋅a_1)⋅a_2 = b_2``.
This gives us ``s = (b_2 - x_1⋅a_2) / (-k a_1⋅a_1)``.
The gap exists if and only if ``s`` is positive.

If the normal vectors do not point in opposite directions, then the defining
hyperplanes intersect and we can produce a witness as follows.
All points ``x`` in this intersection satisfy ``a_1⋅x = b_1`` and
``a_2⋅x = b_2``. Thus we have ``(a_1 + a_2)⋅x = b_1+b_2``.
We now find a dimension where ``a_1 + a_2`` is non-zero, say, ``i``.
Then the result is a vector with one non-zero entry in dimension ``i``, defined
as ``[0, …, 0, (b_1 + b_2)/(a_1[i] + a_2[i]), 0, …, 0]``.
Such a dimension ``i`` always exists.
"""
function isdisjoint(H1::HalfSpace, H2::HalfSpace, witness::Bool=false)
    require(@__MODULE__, :LazySets; fun_name="isdisjoint")

    a1 = H1.a
    a2 = H2.a
    N = promote_type(eltype(H1), eltype(H2))
    issamedir, k = samedir(a1, -a2)
    if issamedir
        x1 = an_element(Hyperplane(a1, H1.b))
        b2 = H2.b
        s = (b2 - dot(x1, a2)) / (-k * dot(a1, a1))
        empty_intersection = s > 0
        # if `!empty_intersection`, x1 is a witness because both defining
        # hyperplanes are contained in each half-space
        return _witness_result_empty(witness, empty_intersection, H1, H2, x1)
    elseif !witness
        return false
    end
    # compute witness
    v = zeros(N, length(a1))
    for i in eachindex(a1)
        a_sum_i = a1[i] + a2[i]
        if !iszero(a_sum_i)
            v[i] = (H1.b + H2.b) / a_sum_i
            break
        end
    end
    return (false, v)
end
