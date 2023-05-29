export box_approximation,
       box_approximation_symmetric,
       symmetric_interval_hull

"""
    box_approximation(S::LazySet{N}) where {N}

Overapproximation a set by a tight hyperrectangle using a parallel
implementation.

### Input

- `S` -- set

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet{N}) where {N}
    (c, r) = box_approximation_helper_parallel(S)
    if r[1] < 0
        return EmptySet{N}(dim(S))
    end
    return Hyperrectangle(c, r)
end

"""
    box_approximation_symmetric(S::LazySet{N}) where {N}

Overapproximate a set by a tight hyperrectangle centered in the origin,
using a parallel implementation.

### Input

- `S` -- set

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N}) where {N}
    (c, r) = box_approximation_helper_parallel(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    box_approximation_helper_parallel(S::LazySet{N}) where {N}

Parallel implementation for the common code of `box_approximation` and
`box_approximation_symmetric`.

### Input

- `S` -- set

### Output

A tuple containing the data that is needed to construct a tightly
overapproximating hyperrectangle.

- `c` -- center
- `r` -- radius

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.

The same load is distributed among all available workers, see
[`distribute_task!`](@ref).
"""
@inline function box_approximation_helper_parallel(S::LazySet{N}) where {N}
    n = dim(S)
    c = SharedVector{N}(n)
    r = SharedVector{N}(n)

    distribute_task!(S, c, r)
    return convert(Vector{N}, c), convert(Vector{N}, r)
end

"""
    process_chunk!(S::LazySet{N}, irange::UnitRange{Int},
                   c::SharedVector{N}, r::SharedVector{N}) where {N}

Kernel to process a given chunk.

### Input

- `c`      -- shared vector with the center of the hyperrectangle
- `r`      -- shared vector with the center of the hyperrectangle
- `S`      -- set
- `irange` -- indices range of the given worker

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.

The load for each worker is passed through the `irange` argument. By default,
the same load is distributed among all available workers. For details see
`distribute_task!`.
"""
function process_chunk!(S::LazySet{N},
                        irange::UnitRange{Int},
                        c::SharedVector{N}, r::SharedVector{N}) where {N}
    d = zeros(N, dim(S))

    for i in irange
        d[i] = one(N)
        htop = ρ(d, S)
        d[i] = -one(N)
        hbottom = -ρ(d, S)
        d[i] = zero(N)
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
    end
end
