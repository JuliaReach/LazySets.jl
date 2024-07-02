"""
    ⊆(B1::Ball2, B2::Ball2, [witness]::Bool=false)

Check whether a ball in the 2-norm is contained in another ball in the 2-norm,
and if not, optionally compute a witness.

### Input

- `B1` -- inner ball in the 2-norm
- `B2` -- outer ball in the 2-norm
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``B1 ⊆ B2``
* If `witness` option is activated:
  * `(true, [])` iff ``B1 ⊆ B2``
  * `(false, v)` iff ``B1 ⊈ B2`` and ``v ∈ B1 ∖ B2``

### Algorithm

``B1 ⊆ B2`` iff ``‖ c_1 - c_2 ‖_2 + r_1 ≤ r_2``
"""
function ⊆(B1::Ball2, B2::Ball2, witness::Bool=false)
    result = norm(B1.center - B2.center, 2) + B1.radius <= B2.radius
    if result
        return _witness_result_empty(witness, true, B1, B2)
    elseif !witness
        return false
    end

    # compute a witness v
    v = B1.center .+ B1.radius * (B1.center .- B2.center)
    return (false, v)
end
