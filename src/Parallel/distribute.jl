export assign_chunk!,
       distribute_task!,
       prange

"""
    assign_chunk!(c::SharedVector{N}, r::SharedVector{N},
                  S::LazySet{N})

Return the function that assigns the work for each process in the
overapproximation of a set by a hyperrectangle.

### Input

- `c` -- center of the hyperrectangle
- `c` -- radius of the hyperrectangle
- `S` -- convex set

### Output

The function `process_chunk!` that equally distributes the load for each worker.
"""
assign_chunk!(c::SharedVector{N},
              r::SharedVector{N},
              S::LazySet{N}) where {N<:Real} = process_chunk!(c, r, S, prange(c))

"""
    distribute_task!(c::SharedVector{N}, r::SharedVector{N},
                     S::LazySet{N}) where {N<:Real}

Distribute the assignment of each chunk among the available processes. 

### Input

- `c` -- center of the hyperrectangle
- `c` -- radius of the hyperrectangle
- `S` -- convex set
"""
function distribute_task!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}) where {N<:Real}
    @sync begin
        for p in procs(c)
            @async remotecall_wait(assign_chunk!, p, c, r, S)
        end
    end
end

"""
    prange(c::SharedVector{N}) where {N<:Real}

Returns the indexes assigned to this process, given a shared vector of
length `n`. 

### Input

- `c` -- shared vector of length `n`

### Output

The tuple `(irange, jrange)` assigned to each process.
"""
function prange(c::SharedVector{N}) where {N<:Real}
    idx = indexpids(c)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(c))
    splits = [round(Int, s) for s in range(0, stop=length(c), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end
