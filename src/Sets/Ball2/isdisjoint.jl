"""
# Extended help

    isdisjoint(B1::Ball2, B2::Ball2, [witness]::Bool=false)

### Algorithm

``B1 ∩ B2 = ∅`` iff ``‖ c_2 - c_1 ‖_2 > r_1 + r_2``.

A witness is computed depending on the smaller/bigger ball (to break ties,
choose `B1` for the smaller ball) as follows.
- If the smaller ball's center is contained in the bigger ball, we return it.
- Otherwise start in the smaller ball's center and move toward the other center
  until hitting the smaller ball's border.
  In other words, the witness is the point in the smaller ball that is closest
  to the center of the bigger ball.
"""
function isdisjoint(B1::Ball2, B2::Ball2, witness::Bool=false)
    center_diff_normed = norm(center(B2) - center(B1), 2)
    empty_intersection = center_diff_normed > B1.radius + B2.radius

    if empty_intersection
        return _witness_result_empty(witness, true, B1, B2)
    elseif !witness
        return false
    end

    # compute a witness 'v' in the intersection
    if B1.radius <= B2.radius
        smaller = B1
        bigger = B2
    else
        smaller = B2
        bigger = B1
    end
    if center_diff_normed <= bigger.radius
        # smaller ball's center is contained in bigger ball
        v = smaller.center
    else
        # scale center difference with smaller ball's radius
        direction = (bigger.center - smaller.center)
        v = smaller.center + direction / center_diff_normed * smaller.radius
    end
    return (false, v)
end
