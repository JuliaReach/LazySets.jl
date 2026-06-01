"""
    constraints_list(B::Bloating)

Return the list of constraints of a bloated set.

### Input

- `B` -- bloated set

### Output

The list of constraints of the bloated set.

### Notes

The constraints list is only available for non-negative bloating in the `p`-norm
for ``p = 1`` or ``p = ∞`` and if `constraints_list` is available for the
unbloated set.

### Algorithm

We call `constraints_list` on the lazy Minkowski sum with the bloating ball.
"""
function constraints_list(B::Bloating)
    @assert ispolyhedral(B) "the constraints list is only available for " *
                            "polyhedral bloating (which requires a polyhedral base set and the " *
                            "1-norm or the infinity norm)"
    if B.ε < 0
        throw(ArgumentError("computing the constraints list of a negatively " *
                            "bloated set is not supported"))
    end

    return constraints_list(MinkowskiSum(B.X, _bloating_ball(B)))
end
