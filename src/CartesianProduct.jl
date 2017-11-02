import Base.*

export CartesianProduct, CartesianProductArray, is_contained

"""
    CartesianProduct <: LazySet

Type that represents the cartesian product.

### Fields

- `X` -- convex set
- `Y` -- another convex set

For the cartesian product a several sets, there exists a special
type `CartesianProductArray`. 
"""
struct CartesianProduct{T1<:LazySet,T2<:LazySet} <: LazySet
    X::T1
    Y::T2
    CartesianProduct{T1,T2}(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = new(X, Y)
    CartesianProduct{T}(Xarr::Vector{T}) where {T<:LazySet} = length(Xarr) == 0 ?
            VoidSet(1) : (length(Xarr) == 1 ? Xarr[1] :
            new{T,T}(Xarr[1], CartesianProduct{T}(Xarr[2:length(Xarr)])))
            # NOTE: use array type instead of element type (bit of a mess otherwise)
end
CartesianProduct(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = CartesianProduct{T1,T2}(X, Y)
CartesianProduct(Xarr::Vector{T}) where {T<:LazySet} = CartesianProduct{T}(Xarr)

function *(X::LazySet, Y::LazySet)::CartesianProduct
    CartesianProduct(X, Y)
end

"""
    dim(cp)

Ambient dimension of a Cartesian product.

### Input

- `cp` -- cartesian product
"""
function dim(cp::CartesianProduct)::Int64
    return dim(cp.X) + dim(cp.Y)
end

"""
    σ(d, cp)

Support vector of a Cartesian product.

### Input

- `d` -- direction

- `cp` -- cartesian product
"""
function σ(d::Vector{Float64}, cp::CartesianProduct)::Vector{Float64}
    return [σ(d[1:dim(cp.X)], cp.X); σ(d[dim(cp.X)+1:end], cp.Y)]
end


"""
    is_contained(d, cp)

Return whether a vector belongs to a given cartesian product set.

### Input

- `d`    --  a vector
- `cp`   -- a cartesian product

### Output

Return true iff d ∈ cp.
"""
function is_contained(d::Vector{Float64}, cp::CartesianProduct)::Bool
    return is_contained(d[1:dim(cp.X)], cp.X) && is_contained(d[dim(cp.X)+1:end], cp.Y)
end

# ==========================
#  Cartesian product of sets
# ==========================
"""
    CartesianProductArray <: LazySet

Type that represents the cartesian product of a finite number of sets.

### Fields

- `sfarray` -- array of sets
"""
mutable struct CartesianProductArray{T<:LazySet} <: LazySet
    sfarray::Vector{T}

    CartesianProductArray{T}(sfarray::Vector{T}) where {T<:LazySet} = new(sfarray)
end
CartesianProductArray() = CartesianProductArray{LazySet}(Vector{LazySet}(0))
CartesianProductArray(sfarray::Vector{T}) where {T<:LazySet} = CartesianProductArray{T}(sfarray)

"""
    dim(cp)

Ambient dimension of the Cartesian product of a finite number of sets.

### Input

- `cp` -- cartesian product array
"""
function dim(cp::CartesianProductArray)::Int64
    return length(cp.sfarray) == 0 ? 0 : sum([dim(sj) for sj in cp.sfarray])
end

"""
    σ(d, cp)

Support vector of the Cartesian product of a finite number of sets.

### Input

- `d` -- direction

- `cp` -- cartesian product array
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

### Input

- `d` -- vector

- `cp` -- cartesian product array
"""
function is_contained(d::Vector{Float64}, cp::CartesianProductArray)::Bool
    contained = false
    jinit = 1
    for Xj in cp
        jend = dim(cp.Xj)
        contained = is_contained(d[jinit:jend], cp.Xj)
        if !contained
            break
        end
        jinit = jend + 1
    end
    return contained
end
