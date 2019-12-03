export volume

"""
    volume(H::AbstractHyperrectangle)

Return the volume of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The (Euclidean) volume of ``H``.

### Algorithm

The volume of the hyperrectangle ``H`` with vector radius ``r`` is
``2ⁿ ∏ᵢ rᵢ`` where ``rᵢ`` denotes the ``i``-th component of ``r``.
"""
function volume(H::AbstractHyperrectangle{N}) where {N<:AbstractFloat}
    r = radius_hyperrectangle(H)
    n = dim(H)
    α = exp(n*log(2))
    vol = α * prod(r)
    return vol
end

# fallback for rational
function volume(H::AbstractHyperrectangle{N}) where {N}
    r = radius_hyperrectangle(H)
    vol = prod(2 .* r)
    return vol
end
