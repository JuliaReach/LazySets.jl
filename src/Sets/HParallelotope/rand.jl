"""
# Extended help

    rand(::Type{HParallelotope}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The directions matrix and offset vector are created randomly. On average there
is a good chance that this resulting set is empty. We then modify the offset to
ensure non-emptiness.

There is a chance that the resulting set represents an unbounded set. This
implementation checks for that case and then samples a new set.
"""
function rand(::Type{HParallelotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    require(@__MODULE__, :LazySets; fun_name="rand")

    rng = reseed!(rng, seed)

    while true
        D = randn(rng, N, dim, dim)
        offset = randn(rng, N, 2 * dim)

        # make sure that the set is not empty:
        # `offset[i] >= -offset[i-dim]` for all `i ∈ dim+1:2*dim`
        @inbounds for i in (dim + 1):(2 * dim)
            offset[i] = -offset[i - dim] + abs(offset[i])
        end

        # convert to polyhedron to check boundedness
        clist = _constraints_list_hparallelotope(D, offset, N, typeof(offset))
        Q = HPolyhedron(clist)
        if isbounded(Q)
            return P = HParallelotope(D, offset; check_consistency=false)
        end
        # set is unbounded; sample a new set in the next iteration
    end  # COV_EXCL_LINE
end
