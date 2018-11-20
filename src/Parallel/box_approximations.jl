export box_approximation,
       box_approximation_symmetric,
       symmetric_interval_hull

"""
    box_approximation(S::LazySet)::Hyperrectangle

Overapproximation a convex set by a tight hyperrectangle using a parallel
algorithm.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet{N};
                          )::Union{Hyperrectangle{N}, EmptySet{N}} where N<:Real
    (c, r) = box_approximation_helper_parallel(S)
    if r[1] < 0
        return EmptySet{N}()
    end
    return Hyperrectangle(c, r)
end

"""
    box_approximation_symmetric(S::LazySet{N}
                               )::Union{Hyperrectangle{N}, EmptySet{N}}
                                where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle centered in the origin,
using a parallel algorithm.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N}) where {N<:Real}
    (c, r) = box_approximation_helper_parallel(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    box_approximation_helper_parallel(S::LazySet{N}) where {N<:Real}

Parallel implementation for the common code of `box_approximation` and
`box_approximation_symmetric`.

### Input

- `S` -- convex set

### Output

A tuple containing the data that is needed to construct a tightly
overapproximating hyperrectangle.

- `c` -- center
- `r` -- radius

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given convex set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.

The same load is distributed among all available workers, see
[`distribute_task!`](@ref).
"""
@inline function box_approximation_helper_parallel(S::LazySet{N}) where {N<:Real}
    n = dim(S)
    c = SharedVector{N}(n)
    r = SharedVector{N}(n)

    distribute_task!(c, r, S)
    return convert(Array, c), convert(Array, r)
end

"""
    process_chunk!(c::SharedVector{N}, r::SharedVector{N},
                   S::LazySet{N}, irange::UnitRange{Int64}) where {N<:Real}

Kernel to process a given chunk 

### Input

- `c`      -- shared vector with the center of the hyperrectangle
- `r`      -- shared vector with the center of the hyperrectangle
- `S`      -- set
- `irange` -- indices range of the given worker

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given convex set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.

The load for each worker is passed through the `irange` argument. By default,
the same load is distributed among all available workers. For details
see `distribute_task!`.
"""
function process_chunk!(c::SharedVector{N}, r::SharedVector{N},
                        S::LazySet{N}, irange::UnitRange{Int64}) where {N<:Real}

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

#=======
Aliases
=======#

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    symmetric_interval_hull

Alias for `box_approximation_symmetric`.
"""
symmetric_interval_hull = box_approximation_symmetric