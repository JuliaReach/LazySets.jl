"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                    ::Type{<:CartesianProductArray{N, S}}
                   ) where {N, S<:LazySet}

Decompose a lazy linear map of a Cartesian product array while keeping the
original block structure.

### Input

- `lm`                    -- lazy linear map of Cartesian product array
- `CartesianProductArray` -- target set type

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N,<:CartesianProductArray},
                         ::Type{<:CartesianProductArray{N,S}}) where {N,S<:LazySet}
    cpa = array(lm.X)
    arr = Vector{S}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, S)
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                    ::Type{<:CartesianProductArray},
                    dir::Type{<:AbstractDirections}) where {N}

Decompose a lazy linear map of a Cartesian product array with template
directions while keeping the original block structure.

### Input

- `lm`                    -- lazy linear map of a Cartesian product array
- `CartesianProductArray` -- target set type
- `dir`                   -- template directions for overapproximation

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N,<:CartesianProductArray},
                         ::Type{<:CartesianProductArray},
                         dir::Type{<:AbstractDirections}) where {N}
    cpa = array(lm.X)
    arr = Vector{HPolytope{N}}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, dir)
end

"""
    overapproximate(lm::LinearMap{N, <:CartesianProductArray},
                    ::Type{<:CartesianProductArray},
                    set_type::Type{<:LazySet}) where {N}

Decompose a lazy linear map of a Cartesian product array with a given set type
while keeping the original block structure.

### Input

- `lm`                    -- lazy linear map of a Cartesian product array
- `CartesianProductArray` -- target set type
- `set_type`              -- set type for overapproximation

### Output

A `CartesianProductArray` representing the decomposed linear map.
"""
function overapproximate(lm::LinearMap{N,<:CartesianProductArray},
                         ::Type{<:CartesianProductArray},
                         set_type::Type{<:LazySet}) where {N}
    cpa = array(lm.X)
    arr = Vector{set_type{N}}(undef, length(cpa))
    return _overapproximate_lm_cpa!(arr, lm.M, cpa, set_type)
end

function _overapproximate_lm_cpa!(arr, M, cpa, overapprox_option)
    # construct Minkowski sum for block row
    function _block_row(cpa::Vector{S}, M::AbstractMatrix,
                        row_range::UnitRange{Int}) where {N,S<:LazySet{N}}
        arr_inner = Vector{LinearMap{N,<:S}}(undef, length(cpa))
        col_start_ind, col_end_ind = 1, 0
        @inbounds for (j, bj) in enumerate(cpa)
            col_end_ind += dim(bj)
            arr_inner[j] = LinearMap(M[row_range, col_start_ind:col_end_ind], bj)
            col_start_ind = col_end_ind + 1
        end
        return MinkowskiSumArray(arr_inner)
    end

    row_start_ind, row_end_ind = 1, 0
    @inbounds for (i, bi) in enumerate(cpa)
        row_end_ind += dim(bi)
        ms = _block_row(cpa, M, row_start_ind:row_end_ind)
        arr[i] = overapproximate(ms, overapprox_option)
        row_start_ind = row_end_ind + 1
    end

    return CartesianProductArray(arr)
end

"""
    overapproximate(rm::ResetMap{N, <:CartesianProductArray},
                    ::Type{<:CartesianProductArray}, oa) where {N}

Overapproximate a reset map (that only resets to zero) of a Cartesian product
with a new Cartesian product.

### Input

- `rm`                    -- reset map
- `CartesianProductArray` -- target set type
- `oa`                    -- overapproximation option

### Output

A Cartesian product with the same block structure.

### Notes

This implementation currently only supports resets to zero.

### Algorithm

We convert the `ResetMap` into a `LinearMap` and then call the corresponding
`overapproximate` method.
"""
function overapproximate(rm::ResetMap{N,<:CartesianProductArray},
                         ::Type{<:CartesianProductArray}, oa) where {N}
    if any(!iszero, values(rm.resets))
        error("this implementation only support resets to zero")
    end

    lm = matrix(rm) * rm.X
    return overapproximate(lm, CartesianProductArray, oa)
end

"""
    overapproximate(cap::Intersection{N,
                                      <:CartesianProductArray,
                                      <:AbstractPolyhedron},
                    ::Type{<:CartesianProductArray}, oa) where {N}

Overapproximate the intersection of a Cartesian product of a finite number of
sets and a polyhedron with a new Cartesian product.

### Input

- `cap`                   -- lazy intersection of a Cartesian product array and
                             a polyhedron
- `CartesianProductArray` -- target set type
- `oa`                    -- overapproximation option

### Output

A `CartesianProductArray` that overapproximates the intersection of `cpa` and
`P`.

### Algorithm

The intersection only needs to be computed in the blocks of `cpa` that are
constrained in `P`.
Hence we first collect those constrained blocks in a lower-dimensional Cartesian
product array and then convert to an `HPolytope` `X`.
Then we take the intersection of `X` and the projection of `Y` onto the
corresponding dimensions.
(This projection is purely syntactic and exact.)
Finally we decompose the result again and plug together the unaffected old
blocks and the newly computed blocks.
The result is a `CartesianProductArray` with the same block structure as in `X`.
"""
function overapproximate(cap::Intersection{N,
                                           <:CartesianProductArray,
                                           <:AbstractPolyhedron},
                         ::Type{<:CartesianProductArray}, oa) where {N}
    cpa, P = first(cap), second(cap)

    cpa_low_dim, vars, block_structure, blocks = get_constrained_lowdimset(cpa, P)

    hpoly_low_dim = HPolytope(constraints_list(cpa_low_dim))
    low_intersection = intersection(hpoly_low_dim, project(P, vars))

    if isempty(low_intersection)
        return EmptySet{N}(dim(cap))
    end

    decomposed_low_set = decompose(low_intersection, block_structure, oa)

    return substitute_blocks(decomposed_low_set, cpa, blocks)
end

# symmetric method
function overapproximate(cap::Intersection{N,
                                           <:AbstractPolyhedron,
                                           <:CartesianProductArray},
                         ::Type{<:CartesianProductArray}, oa) where {N}
    return overapproximate(swap(cap), oa)
end
