@validate function project(cpa::CartesianProductArray, block::AbstractVector{Int}; kwargs...)
    target_sets = LazySet[]
    m = length(block)

    # find first set
    i_start = 1
    bi = block[i_start]
    n_sum = 0
    n_sum_old = 0
    @inbounds for (j, Xj) in enumerate(array(cpa))
        nj = dim(Xj)
        n_sum += nj
        if n_sum >= bi
            # found starting point in a set; now find end point
            i_end = m
            for i in (i_start + 1):m
                if block[i] > n_sum
                    i_end = i - 1
                    break
                end
            end

            # project this block
            projected = project(Xj, block[i_start:i_end] .- n_sum_old; kwargs...)
            push!(target_sets, projected)

            if i_end == m
                # last index visited
                break
            end

            # advance indices
            i_start = i_end + 1
            bi = block[i_start]
        end
        n_sum_old = n_sum
    end

    # construct result depending on the number of sets
    if length(target_sets) == 1
        @inbounds return target_sets[1]
    elseif length(target_sets) == 2
        @inbounds return CartesianProduct(target_sets[1], target_sets[2])
    else
        # create a new array for better type information
        return CartesianProductArray([X for X in target_sets])
    end
end
