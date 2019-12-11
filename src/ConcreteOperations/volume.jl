export volume

"""
    volume(H::AbstractHyperrectangle{N}) where {N<:Real}

Return the volume of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The volume of ``H``.

### Algorithm

The volume of the ``n``-dimensional hyperrectangle ``H`` with vector radius ``r`` is
``2ⁿ ∏ᵢ rᵢ`` where ``rᵢ`` denotes the ``i``-th component of ``r``.
"""
function volume(H::AbstractHyperrectangle{N}) where {N<:Real}
    vol = mapreduce(x -> 2x, *, radius_hyperrectangle(H))
    return vol
end

function volume(B::BallInf{N}) where {N<:Real}
    n = dim(B)
    if n < 50
        vol = one(N)
        for i in 1:n
            vol *= 2 * radius_hyperrectangle(B, i)
        end
    else
        r = B.radius
        vol = exp(n*log(2r))
    end
    return vol
end
