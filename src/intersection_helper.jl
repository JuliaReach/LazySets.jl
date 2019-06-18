"""
    get_constrained_lowdimset(cpa::CartesianProductArray{N, S}, P::AbstractPolyhedron{N}) where {N<:Real, S<:LazySet{N}}

Preprocess step for intersection between cartesian product array and H-polyhedron.
Returns low-dimensional HPolytope in constrained dimensions of original cpa,
constrained variables and variables in corresponding blocks, original block structure
of low-dimensional set and list of constrained blocks.

### Input

- `cpa` -- Cartesian product array of convex sets
- `P` -- polyhedron

### Output

A tuple of low-dimensional set, list of constrained dimensions, original block
structure of low-dimensional set and corresponding blocks indices.
"""
function get_constrained_lowdimset(cpa::CartesianProductArray{N, S}, P::AbstractPolyhedron{N}) where {N<:Real, S<:LazySet{N}}

    if isbounded(P)
        blocks, non_empty_length = block_to_dimension_indices(cpa)
    else
        constrained_vars = constrained_dimensions(P)
        blocks, non_empty_length = block_to_dimension_indices(cpa, constrained_vars)
    end

    array = Vector{S}()
    sizehint!(array, non_empty_length)
    vars = Vector{Int}()
    block_structure = Vector{UnitRange{Int}}()
    sizehint!(block_structure, non_empty_length)

    last_var = 1
    for i in 1:length(blocks)
        start_index, end_index = blocks[i]
        block_end = last_var + end_index - start_index
        if start_index != -1
            push!(array, cpa.array[i])
            append!(vars, start_index : end_index)
            push!(block_structure, last_var : block_end)
            last_var = block_end + 1
        end
    end

    cpa_low_dim = HPolytope(constraints_list(CartesianProductArray(array)));

    return cpa_low_dim, vars, block_structure, blocks
end
