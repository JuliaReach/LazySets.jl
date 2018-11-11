export assign_chunk!,
       distribute_task!,
       myrange

assign_chunk!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}) where {N<:Real} = process_chunk!(c, r, S, myrange(c))

function distribute_task!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}) where {N<:Real}
    @sync begin
        for p in procs(c)
            @async remotecall_wait(assign_chunk!, p, c, r, S)
        end
    end
end

# This function retuns the (irange,jrange) indexes assigned to this worker
function myrange(c::SharedVector{N}) where {N<:Real}
    idx = indexpids(c)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(c))
    splits = [round(Int, s) for s in range(0, stop=length(c), length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end
