"""
# Extended help

    ρ(d::AbstractVector, B::Ball1)

### Algorithm

Let ``c`` and ``r`` be the center and radius of the ball ``B`` in the 1-norm,
respectively. Then:

```math
ρ(d, B) = ⟨d, c⟩ + r ‖d‖_∞.
```
"""
function ρ(d::AbstractVector, B::Ball1)
    return dot(d, B.center) + B.radius * maximum(abs, d)
end
