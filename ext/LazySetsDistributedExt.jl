module LazySetsDistributedExt

using Distributed: remotecall_wait, procs
using LazySets: LazySet, dim, ρ
using LazySets.EmptySetModule: EmptySet
using LazySets.HyperrectangleModule: Hyperrectangle
using SharedArrays: SharedVector, indexpids

export assign_chunk!, distribute_task!,
       box_approximation, box_approximation_symmetric, symmetric_interval_hull

#=======================================================
Utility functions for distribution of tasks in parallel
=======================================================#

"""
    distribute_task!(X::LazySet, v::SharedVector...)

Distribute the assignment of each chunk among the available processes.

### Input

- `X` -- set
- `v` -- shared vectors

### Output

Nothing.

### Notes

Use this function to distribute a given task acting on a set `X` and a pool `v`
of shared vectors.

The task for each processor is distributed through `remotecall_wait` using a
function `assign_chunk!` that should be defined elsewhere. The vectors `v`
contain one or more shared vectors in which the values of the task are written.
"""
function distribute_task!(X::LazySet, v::SharedVector...)
    @sync begin
        for p in procs(v[1])
            @async remotecall_wait(assign_chunk!, p, X, v...)
        end
    end
end

"""
    assign_chunk!(X::LazySet, v::SharedVector...)

Return the function that assigns the work for each process.

### Input

- `X` -- set
- `v` -- shared vectors

### Output

The function `process_chunk!` that equally distributes the load for each worker.

### Notes

This function is a wrapper around a problem-specific `process_chunk!` function.

Use this function to distribute a given task acting on a set `X` and a pool `v`
of shared vectors. The tasks are equally distributed among the processes.

See also [`distribute_task!`](@ref).
"""
function assign_chunk!(X::LazySet, v::SharedVector...)
    return process_chunk!(X, _prange(v[1]), v...)
end

"""
    _prange(v::SharedVector)

Returns the indices assigned to a process.

### Input

- `v` -- shared vector

### Output

The indices assigned to each process.

### Notes

The indices are assigned such that the vector is equally distributed among the
processes. If the worker is not assigned a piece, the unit range `1:0` is
returned.
"""
function _prange(v::SharedVector)
    idx = indexpids(v)
    if idx == 0
        # this worker is not assigned a piece
        return 1:0
    end
    nchunks = length(procs(v))
    splits = [round(Int, s) for s in range(0; stop=length(v), length=nchunks + 1)]
    return (splits[idx] + 1):splits[idx + 1]
end

#==================================================
Approximations using boxes implemented in parallel
==================================================#

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

end  # module
