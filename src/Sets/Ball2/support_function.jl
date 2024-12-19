"""
# Extended help

    ρ(d::AbstractVector, B::Ball2)

### Algorithm

Let ``c`` and ``r`` be the center and radius of the ball ``B`` in the 2-norm,
respectively. Then:

```math
ρ(d, B) = ⟨d, c⟩ + r ‖d‖_2.
```
"""
function ρ(d::AbstractVector, B::Ball2)
    return dot(d, B.center) + B.radius * norm(d, 2)
end
