const BALLINF_THRESHOLD_VOLUME = 50  # threshold value in `volume`

"""
    volume(B::BallInf)

Return the volume of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The volume of ``B``.

### Algorithm

We compute the volume by iterative multiplication of the radius.

For floating-point inputs we use this implementation for balls of dimension less
than ``$BALLINF_THRESHOLD_VOLUME``. For balls of higher dimension we instead compute
``exp(n * log(2r))``, where ``r`` is the radius of the ball.
"""
function volume(B::BallInf)
    return _pow_loop(B.radius, dim(B))
end

# method for floating-point input
function volume(B::BallInf{N}) where {N<:AbstractFloat}
    n = dim(B)
    if n < BALLINF_THRESHOLD_VOLUME
        vol = _pow_loop(B.radius, n)
    else
        vol = _pow_exp(B.radius, n)
    end
    return vol
end

# compute a^n in a loop
@inline function _pow_loop(a::N, n::Int) where {N}
    vol = one(N)
    diam = 2 * a
    for i in 1:n
        vol *= diam
    end
    return vol
end

# compute a^n using exp
@inline function _pow_exp(a, n::Int)
    return exp(n * log(2a))
end
