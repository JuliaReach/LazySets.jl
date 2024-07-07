"""
    ρ(d::AbstractVector, B::Ball1)

Evaluate the support function of a ball in the 1-norm in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the 1-norm

### Output

Evaluation of the support function in the given direction.

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
