# the algorithm is a 2D specialization of the `Hyperplane` algorithm, except
# that it returns a `Singleton` for a 1D line
function project(L::Line2D{N}, block::AbstractVector{Int}; kwargs...) where {N}
    m = length(block)
    if m == 2
        @inbounds if block[1] == 1 && block[2] == 2
            return L  # no projection
        elseif block[1] == 2 && block[2] == 1
            return Line2D(L.a[block], L.b)  # swap a vector
        else
            throw(ArgumentError("invalid projection to $block"))
        end
    elseif m == 1
        # projection to dimension i
        cdims = constrained_dimensions(L)
        if length(cdims) == 1
            @inbounds if cdims[1] == block[1]
                # L: aᵢxᵢ = b where aᵢ ≠ 0
                return Singleton([L.b / L.a[cdims[1]]])
            else
                # L: aⱼxⱼ = b where i ≠ j
                return Universe{N}(1)
            end
        else
            # L is constrained in both dimensions
            @assert length(cdims) == 2
            return Universe{N}(1)
        end
    else
        throw(ArgumentError("cannot project a two-dimensional line to $m dimensions"))
    end
end

"""
    project(x::AbstractVector, L::Line2D)

Project a point onto a 2D line.

### Input

- `x` -- point/vector
- `L` -- 2D line

### Output

The projection of `x` onto `L`.

### Algorithm

The projection of ``x`` onto a line of the form ``a⋅x = b`` is

```math
    x - \\dfrac{a (a⋅x - b)}{‖a‖²}.
```
"""
function project(x::AbstractVector, L::Line2D)
    return x - L.a * (dot(L.a, x) - L.b) / norm(L.a, 2)^2
end
