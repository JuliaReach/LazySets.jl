export RotatedHyperrectangle

"""
    RotatedHyperrectangle{N, MN<:AbstractMatrix{N},
                          HT<:AbstractHyperrectangle{N}} <: AbstractZonotope{N}

Type that represents a hyperrectangle that is not necessarily axis aligned.

### Fields

- `M`   -- matrix (not necessarily invertible)
- `box` -- axis-aligned hyperrectangle

### Notes

The matrix `M` is typically a rotation matrix and hence invertible, but this
type does not require that and any matrix that is compatible with the dimension
of `box` is allowed.
"""
struct RotatedHyperrectangle{N, MN<:AbstractMatrix{N},
                             HT<:AbstractHyperrectangle{N}} <: AbstractZonotope{N}
    M::MN
    box::HT

    function RotatedHyperrectangle(M::MN, box::HT) where {N,
            MN<:AbstractMatrix{N}, HT<:AbstractHyperrectangle{N}}
        @assert size(M, 2) == dim(box) "a hyperrectangle of dimension " *
            "$(dim(box)) is incompatible with a matrix of dimension $(size(M))"
        return new{N, MN, HT}(M, box)
    end
end

isoperationtype(::Type{RotatedHyperrectangle}) = false

"""
    dim(R::RotatedHyperrectangle)

Return the ambient dimension of a rotated hyperrectangle.

### Input

- `R` -- rotated hyperrectangle

### Output

The ambient dimension of the rotated hyperrectangle.
"""
function dim(R::RotatedHyperrectangle)
    return size(R.M, 1)
end

"""
    ρ(d::AbstractVector, R::RotatedHyperrectangle)

Evaluate the support function of a rotated hyperrectangle in a given direction.

### Input

- `d` -- direction
- `R` -- rotated hyperrectangle

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, R::RotatedHyperrectangle)
    return _ρ_linear_map(d, R.M, R.box)
end

"""
    σ(d::AbstractVector, R::RotatedHyperrectangle)

Return a support vector of a rotated hyperrectangle in a given direction.

### Input

- `d` -- direction
- `R` -- rotated hyperrectangle

### Output

A support vector in the given direction.
"""
function σ(d::AbstractVector, R::RotatedHyperrectangle)
    return _σ_linear_map(d, R.M, R.box)
end

"""
    center(R::RotatedHyperrectangle)

Return the center of a rotated hyperrectangle.

### Input

- `R` -- rotated hyperrectangle

### Output

The center of the rotated hyperrectangle.
"""
function center(R::RotatedHyperrectangle)
    return R.M * center(R.box)
end

"""
    generators(R::RotatedHyperrectangle)

Return an iterator over the generators of a rotated hyperrectangle.

### Input

- `R` -- rotated hyperrectangle

### Output

An iterator over the generators of `R`.
"""
function generators(R::RotatedHyperrectangle)
    return generators_fallback(R)
end

"""
   genmat(R::RotatedHyperrectangle)

Return the generator matrix of a rotated hyperrectangle.

### Input

- `R` -- rotated hyperrectangle

### Output

A matrix where each column represents one generator of `R`.
"""
function genmat(R::RotatedHyperrectangle)
    return R.M * genmat(R.box)
end

"""
    linear_map(M::AbstractMatrix, R::RotatedHyperrectangle)

Compute the concrete linear map of a rotated hyperrectangle.

### Input

- `M` -- matrix
- `R` -- rotated hyperrectangle

### Output

A new rotated hyperrectangle.

### Notes

If `M` is not a square matrix, the result is not necessarily a rotated
hyperrectangle.
The type still represents the correct set, but it is suggested to only apply
square linear maps to it.
"""
function linear_map(M::AbstractMatrix, R::RotatedHyperrectangle)
    @assert size(M, 2) == dim(R) "a linear map of size $(size(M)) is " *
        "incompatible with a set of dimension $(dim(R))"

    return RotatedHyperrectangle(M * R.M, R.box)
end

"""
    vertices_list(R::RotatedHyperrectangle)

Return the list of vertices of a rotated hyperrectangle.

### Input

- `R` -- rotated hyperrectangle

### Output

A list of the vertices.
"""
function vertices_list(R::RotatedHyperrectangle)
    return broadcast(v -> R.M * v, vertices_list(R.box))
end

"""
    constraints_list(R::RotatedHyperrectangle)

Return the list of constraints of a rotated hyperrectangle.

### Input

- `R` -- rotated hyperrectangle

### Output

A list of the linear constraints.
"""
function constraints_list(R::RotatedHyperrectangle)
    return constraints_list(convert(Zonotope, R))
end
