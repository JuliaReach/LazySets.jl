"""
    get_constrained_lowdimset(cpa::CartesianProductArray{N, S},
                              P::AbstractPolyhedron{N}
                             ) where {N<:Real, S<:LazySet{N}}

Preprocess step for intersection between Cartesian product array and polyhedron.
Returns low-dimensional a `CartesianProductArray` in the constrained dimensions
of the original cpa,
constrained variables and variables in corresponding blocks, original block
structure of low-dimensional set and list of constrained blocks.

### Input

- `cpa` -- Cartesian product array of convex sets
- `P`   -- polyhedron

### Output

A tuple of low-dimensional set, list of constrained dimensions, original block
structure of low-dimensional set and corresponding blocks indices.
"""
function get_constrained_lowdimset(cpa::CartesianProductArray{N, S},
                                   P::AbstractPolyhedron{N}
                                  ) where {N<:Real, S<:LazySet{N}}

    if isbounded(P)
        blocks, non_empty_length = block_to_dimension_indices(cpa)
    else
        blocks, non_empty_length =
            block_to_dimension_indices(cpa, constrained_dimensions(P))
    end

    array = Vector{S}()
    sizehint!(array, non_empty_length)
    variables = Vector{Int}()
    block_structure = Vector{UnitRange{Int}}()
    sizehint!(block_structure, non_empty_length)

    last_var = 1
    for i in 1:length(blocks)
        start_index, end_index = blocks[i]
        block_end = last_var + end_index - start_index
        if start_index != -1
            push!(array, cpa.array[i])
            append!(variables, start_index : end_index)
            push!(block_structure, last_var : block_end)
            last_var = block_end + 1
        end
    end

    return CartesianProductArray(array), variables, block_structure, blocks
end
