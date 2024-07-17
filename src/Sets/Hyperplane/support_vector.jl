"""
    σ(d::AbstractVector, H::Hyperplane)

Return a support vector of a hyperplane.

### Input

- `d` -- direction
- `H` -- hyperplane

### Output

A support vector in the given direction, which is only defined in the following
two cases:
1. The direction has norm zero.
2. The direction is the hyperplane's normal direction or its opposite direction.
In all cases, any point on the hyperplane is a solution.
Otherwise this function throws an error.
"""
function σ(d::AbstractVector, H::Hyperplane)
    v, unbounded = _σ_hyperplane_halfspace(d, H.a, H.b; error_unbounded=true,
                                           halfspace=false)
    return v
end

"""
```
    _σ_hyperplane_halfspace(d::AbstractVector, a, b;
                            [error_unbounded]::Bool=true,
                            [halfspace]::Bool=false)
```

Return a support vector of a hyperplane ``a⋅x = b`` in direction `d`.

### Input

- `d`         -- direction
- `a`         -- normal direction
- `b`         -- constraint
- `error_unbounded` -- (optional, default: `true`) `true` if an error should be
                 thrown whenever the set is unbounded in the given direction
- `halfspace` -- (optional, default: `false`) `true` if the support vector
                 should be computed for a half-space

### Output

A pair `(v, f)` where `v` is a vector and `f` is a Boolean flag.

The flag `f` is `false` in one of the following cases:
1. The direction has norm zero.
2. The direction is (a multiple of) the hyperplane's normal direction.
3. The direction is (a multiple of) the opposite of the hyperplane's normal
direction and `halfspace` is `false`.
In all these cases, `v` is any point on the hyperplane.

Otherwise, the flag `f` is `true`, the set is unbounded in the given direction,
and `v` is any vector.

If `error_unbounded` is `true` and the set is unbounded in the given direction,
this function throws an error instead of returning.

### Notes

For correctness, consider the
[weak duality of LPs](https://en.wikipedia.org/wiki/Linear_programming#Duality):
If the primal is unbounded, then the dual is infeasible.
Since there is only a single constraint, the feasible set of the dual problem is
``a ⋅ y == d``, ``y ≥ 0`` (with objective function ``b ⋅ y``).
It is easy to see that this problem is infeasible whenever ``a`` is not parallel
to ``d``.
"""
@inline function _σ_hyperplane_halfspace(d::AbstractVector, a, b;
                                         error_unbounded::Bool=true,
                                         halfspace::Bool=false)
    @assert length(d) == length(a) "cannot compute the support vector of a " *
                                   "$(length(a))-dimensional " *
                                   (halfspace ? "halfspace" : "hyperplane") *
                                   " along a vector of length $(length(d))"

    first_nonzero_entry_a = -1
    unbounded = false
    if iszero(d)
        # zero vector
        return (_an_element_helper_hyperplane(a, b), false)
    else
        # not the zero vector, check if it is a normal vector
        N = promote_type(eltype(d), eltype(a))
        factor = zero(N)
        for i in eachindex(a)
            if a[i] == 0
                if d[i] != 0
                    unbounded = true
                    break
                end
            else
                if d[i] == 0
                    unbounded = true
                    break
                elseif first_nonzero_entry_a == -1
                    factor = a[i] / d[i]
                    first_nonzero_entry_a = i
                    if halfspace && factor < 0
                        unbounded = true
                        break
                    end
                elseif d[i] * factor != a[i]
                    unbounded = true
                    break
                end
            end
        end
        if !unbounded
            return (_an_element_helper_hyperplane(a, b, first_nonzero_entry_a), false)
        end
        if error_unbounded
            error("the support vector for the " *
                  (halfspace ? "halfspace" : "hyperplane") * " with normal " *
                  "direction $a is not defined along a direction $d")
        end
        # the first return value does not have a meaning here
        return (d, true)
    end
end
