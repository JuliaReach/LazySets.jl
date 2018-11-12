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

function box_approximation_symmetric(S::LazySet{N}) where {N<:Real}
    (c, r) = box_approximation_helper_parallel(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

@inline function box_approximation_helper_parallel(S::LazySet{N}) where {N<:Real}
    n = dim(S)
    c = SharedVector{N}(n)
    r = SharedVector{N}(n)

    distribute_task!(c, r, S)
    return convert(Array, c), convert(Array, r)
end

# Here's the kernel
function process_chunk!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}, irange::UnitRange{Int64}) where {N<:Real}

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