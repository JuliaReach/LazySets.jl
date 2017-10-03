"""
    CartesianProduct <: LazySet

Type that represents the cartesian product.

FIELDS:

- ``s1`` -- convex set
- ``s2`` -- another convex set

For the cartesian product a several sets, it is recommended to use the dedicated
type ``CartesianProductArray``. 
"""
struct CartesianProduct <: LazySet
    s1::LazySet
    s2::LazySet
    CartesianProduct(s1::LazySet, s2::LazySet) = new(s1, s2)
    CartesianProduct(as::Array{LazySet, 1}) = length(as) == 0 ? VoidSet(1) : (length(as) == 1 ? as[1] : new(as[1], CartesianProduct(as[2:length(as)])))
end


import Base: *

function *(s1::LazySet, s2::LazySet)::CartesianProduct
    CartesianProduct(s1, s2)
end

"""
    dim(cp)

Ambient dimension of a Cartesian product.

INPUT:

- ``cp`` -- cartesian product
"""
function dim(cp::CartesianProduct)::Int64
    return dim(cp.s1) + dim(cp.s2)
end

"""
    σ(d, cp)

Support vector of a Cartesian product.

INPUT:

- ``d`` -- direction

- ``cp`` -- cartesian product
"""
function σ(d::Vector{Float64}, cp::CartesianProduct)::Vector{Float64}
    return [σ(d[1:dim(cp.s1)], cp.s1); σ(d[dim(cp.s1)+1:end], cp.s2)]
end


"""
    is_contained(d, cp)

Return whether a vector belongs to a given cartesian product set.

INPUT :

- ``d``    --  a vector
- ``cp``   -- a cartesian product

OUTPUT :

Return true iff d ∈ cp.
"""
function is_contained(d::Vector{Float64}, cp::CartesianProduct)::Bool
    return is_contained(d[1:dim(cp.s1)], cp.s1) && is_contained(d[dim(cp.s1)+1:end], cp.s2)
end

export CartesianProduct


# ============ Cartesian product of sets ================
"""
    CartesianProductArray <: LazySet

Type that represents the cartesian product of a finite number of sets.

FIELDS:

- ``sfarray`` -- array of sets
"""
mutable struct CartesianProductArray <: LazySet
    sfarray::Array{LazySet, 1}

    CartesianProductArray(sfarray) = new(sfarray)
end
CartesianProductArray() = CartesianProductArray(Array{LazySet, 1}(0))

"""
    dim(cp)

Ambient dimension of the Cartesian product of a finite number of sets.

INPUT:

- ``cp`` -- cartesian product array
"""
function dim(cp::CartesianProductArray)::Int64
    return length(cp.sfarray) == 0 ? 0 : sum([dim(sj) for sj in cp.sfarray])
end

"""
    σ(d, cp)

Support vector of the Cartesian product of a finite number of sets.

INPUT:

- ``d`` -- direction

- ``cp`` -- cartesian product array
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, cp::CartesianProductArray)::Vector{Float64}
    svec = Vector{Float64}(length(d))
    jinit = 1
    for sj in cp.sfarray
        jend = jinit + dim(sj) - 1
        svec[jinit:jend] = σ(d[jinit:jend], sj)
        jinit = jend + 1
    end
    return svec
end

"""
    is_contained(d, cp)

Return whether a given vector is contained in the cartesian product of a
finite number of sets.

INPUT:

- ``d`` -- vector

- ``cp`` -- cartesian product array
"""
function is_contained(d::Vector{Float64}, cp::CartesianProductArray)::Bool
    contained = false
    jinit = 1
    for sj in cp
        jend = dim(cp.sj)
        contained = is_contained(d[jinit:jend], cp.sj)
        if !contained
            break
        end
        jinit = jend + 1
    end
    return contained
end

export CartesianProductArray, is_contained
