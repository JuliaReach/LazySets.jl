export assign_chunk!,
       distribute_task!

"""
    distribute_task!(S::LazySet{N}, v::SharedVector{N}...) where {N<:Real}

Distribute the assignment of each chunk among the available processes.

### Input

- `S` -- convex set
- `v` -- variable number of shared vectors

### Output

Nothing.

### Notes

Use this function to distribute a given task acting on a set `S` and a pool `v`
of shared vectors.

The task for each processor is distributed through `remotecall_wait` using a
function `assign_chunk!` that should be defined elsewhere. The vectors `v`
contain one or more shared vectors in which the values of the task are written.

Typically, the function `assign_chunk!` is a wrapper around some problem-specific
`process_chunk!` function.
"""
function distribute_task!(S::LazySet{N}, v::SharedVector{N}...) where {N<:Real}
    @sync begin
        for p in procs(v[1])
            @async remotecall_wait(assign_chunk!, p, S, v...)
        end
    end
end

"""
    assign_chunk!(S::LazySet{N}, v::SharedVector{N}...) where {N<:Real}

Return the function that assigns the work for each process.

### Input

- `S` -- convex set
- `v` -- variable number of shared vectors

### Output

The function `process_chunk!` that equally distributes the load for each worker.

### Notes

Use this function to distribute a given task acting on a set `S` and a pool `v`
of shared vectors. The tasks are equally distributed among the number of processes.

See also [`distribute_task!`](@ref).
"""
function assign_chunk!(S::LazySet{N}, v::SharedVector{N}...) where {N<:Real}
    return process_chunk!(S, _prange(v[1]), v...)
end

"""
    _prange(v::SharedVector{N}) where {N<:Real}

Returns the indexes assigned to a process.

### Input

- `v` -- shared vector of length `n`

### Output

The indices range assigned to each process.

### Notes

The indices are assigned such that the vector is equally distributed among the
processes. If the worker is not assigned a piece, the unit range `1:0` is
returned.
"""
function _prange(v::SharedVector{N}) where {N<:Real}
    idx = indexpids(v)
    if idx == 0
        # this worker is not assigned a piece
        return 1:0
    end
    nchunks = length(procs(v))
    splits = [round(Int, s) for s in range(0, stop=length(v), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end
