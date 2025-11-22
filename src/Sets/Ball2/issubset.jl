"""
# Extended help

    issubset(B1::Ball2, B2::Ball2, [witness]::Bool=false)

### Algorithm

``B1 ⊆ B2`` iff ``‖ c_1 - c_2 ‖_2 + r_1 ≤ r_2``
"""
@validate function issubset(B1::Ball2, B2::Ball2, witness::Bool=false)
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
