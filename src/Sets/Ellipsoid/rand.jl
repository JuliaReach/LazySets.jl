"""
    rand(::Type{Ellipsoid}; [N]::Type{<:AbstractFloat}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random ellipsoid.

### Input

- `Ellipsoid` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random ellipsoid.

### Algorithm

The center is a normally distributed vector with entries of mean 0 and standard
deviation 1.

The idea for the shape matrix comes from
[here](https://math.stackexchange.com/a/358092).
The matrix is symmetric positive definite, but also diagonally dominant.

```math
Q =  \\frac{1}{2}(S + S^T) + nI,
```
where ``n`` = `dim` and ``S`` is a ``n × n`` random matrix whose
coefficients are uniformly distributed in the interval ``[-1, 1]``.
"""
function rand(::Type{Ellipsoid};
              N::Type{<:AbstractFloat}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    center = randn(rng, N, dim)
    # random entries in [-1, 1]
    # this needs a bit of code because 'rand' only samples from [0, 1]
    shape_matrix = Matrix{N}(undef, dim, dim)
    for j in 1:dim
        for i in 1:dim
            entry = rand(rng, N)
            if rand(rng, Bool)
                entry = -entry
            end
            shape_matrix[i, j] = entry
        end
    end
    # make diagonally dominant
    shape_matrix = N(0.5) * (shape_matrix + shape_matrix') +
                   Matrix{N}(dim * I, dim, dim)
    return Ellipsoid(center, shape_matrix)
end
