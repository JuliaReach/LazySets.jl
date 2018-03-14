import Base:*, ×, ∈

export CartesianProduct,
       CartesianProductArray,
       CartesianProduct!,
       array

"""
    CartesianProduct{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents a Cartesian product of two convex sets.

### Fields

- `X` -- first convex set
- `Y` -- second convex set

### Notes

The Cartesian product of three elements is obtained recursively.
See also `CartesianProductArray` for an implementation of a Cartesian product of
many sets without recursion, instead using an array.

The `EmptySet` is the absorbing element for `CartesianProduct`.

Constructors:

- `CartesianProduct{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}(X1::S1, X2::S2)`
  -- default constructor

- `CartesianProduct(Xarr::Vector{S}) where {S<:LazySet}`
  -- constructor from an array of convex sets
"""
struct CartesianProduct{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
end
# type-less convenience constructor
CartesianProduct(X1::S1, X2::S2
                ) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    CartesianProduct{N, S1, S2}(X1, X2)
# constructor from an array
CartesianProduct(Xarr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
    (length(Xarr) == 0
        ? EmptySet{N}()
        : length(Xarr) == 1
            ? Xarr[1]
            : length(Xarr) == 2
                ? CartesianProduct(Xarr[1], Xarr[2])
                : CartesianProduct(Xarr[1],
                                   CartesianProduct(Xarr[2:length(Xarr)])))

# EmptySet is the absorbing element for CartesianProduct
@absorbing(CartesianProduct, EmptySet)

"""
```
    *(X::LazySet, Y::LazySet)
```

Alias for the binary Cartesian product.
"""
*(X::LazySet, Y::LazySet) = CartesianProduct(X, Y)

"""
    ×

Alias for the binary Cartesian product.
"""
×(X::LazySet, Y::LazySet) = CartesianProduct(X, Y)

"""
    dim(cp::CartesianProduct)::Int

Return the dimension of a Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

The ambient dimension of the Cartesian product.
"""
function dim(cp::CartesianProduct)::Int
    return dim(cp.X) + dim(cp.Y)
end

"""
    σ(d::V, cp::CartesianProduct) where {N<:Real, V<:AbstractVector{N}}

Return the support vector of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the product sets.

### Algorithm


"""
function σ(d::V, cp::CartesianProduct) where {N<:Real, V<:AbstractVector{N}}
    return [σ(d[1:dim(cp.X)], cp.X); σ(d[dim(cp.X)+1:length(d)], cp.Y)]
end

"""
    ∈(x::AbstractVector{<:Real}, cp::CartesianProduct)::Bool

Check whether a given point is contained in a Cartesian product set.

### Input

- `x`  -- point/vector
- `cp` -- Cartesian product

### Output

`true` iff ``x ∈ cp``.
"""
function ∈(x::AbstractVector{<:Real}, cp::CartesianProduct)::Bool
    @assert length(x) == dim(cp)

    return ∈(view(x, 1:dim(cp.X)), cp.X) &&
           ∈(view(x, dim(cp.X)+1:length(x)), cp.Y)
end

# ======================================
#  Cartesian product of an array of sets
# ======================================
"""
    CartesianProductArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the Cartesian product of a finite number of convex sets.

### Fields

- `array` -- array of sets

### Notes

The `EmptySet` is the absorbing element for `CartesianProductArray`.

Constructors:

- `CartesianProductArray(array::Vector{<:LazySet})` -- default constructor

- `CartesianProductArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty product with optional size hint and numeric type
"""
struct CartesianProductArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

# type-less convenience constructor
CartesianProductArray(arr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
    CartesianProductArray{N, S}(arr)

# constructor for an empty product with optional size hint and numeric type
function CartesianProductArray(n::Int=0, N::Type=Float64)::CartesianProductArray
    arr = Vector{LazySet{N}}(0)
    sizehint!(arr, n)
    return CartesianProductArray(arr)
end

# EmptySet is the absorbing element for CartesianProductArray
@absorbing(CartesianProductArray, EmptySet)

# add functions connecting CartesianProduct and CartesianProductArray
@declare_array_version(CartesianProduct, CartesianProductArray)

"""
    array(cpa::CartesianProductArray{N, S}
         )::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a Cartesian product of a finite number of convex sets.

### Input

- `cpa` -- Cartesian product array

### Output

The array of a Cartesian product of a finite number of convex sets.
"""
function array(cpa::CartesianProductArray{N, S}
              )::Vector{S} where {N<:Real, S<:LazySet{N}}
    return cpa.array
end

"""
    dim(cpa::CartesianProductArray)::Int

Return the dimension of a Cartesian product of a finite number of convex sets.

### Input

- `cpa` -- Cartesian product array

### Output

The ambient dimension of the Cartesian product of a finite number of convex
sets.
"""
function dim(cpa::CartesianProductArray)::Int
    return length(cpa.array) == 0 ? 0 : sum([dim(sj) for sj in cpa.array])
end

"""
    σ(d::V, cpa::CartesianProductArray{N, <:LazySet{N}}) where {N<:Real, V<:AbstractVector{N}}

Support vector of a Cartesian product.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the product sets.
"""
function σ(d::V, cpa::CartesianProductArray{N, <:LazySet{N}}) where {N<:Real, V<:AbstractVector{N}}
    svec = similar(d)
    jinit = 1
    for sj in cpa.array
        jend = jinit + dim(sj) - 1
        svec[jinit:jend] = σ(d[jinit:jend], sj)
        jinit = jend + 1
    end
    return svec
end

"""
    ∈(x::AbstractVector{N}, cpa::CartesianProductArray{N, <:LazySet{N}}
     )::Bool  where {N<:Real}

Check whether a given point is contained in a Cartesian product of a finite
number of sets.

### Input

- `x`   -- point/vector
- `cpa` -- Cartesian product array

### Output

`true` iff ``x ∈ \\text{cpa}``.
"""
function ∈(x::AbstractVector{N}, cpa::CartesianProductArray{N, <:LazySet{N}}
          )::Bool  where {N<:Real}
    @assert length(x) == dim(cpa)

    jinit = 1
    for sj in cpa.array
        jend = jinit + dim(sj) - 1
        if !∈(x[jinit:jend], sj)
            return false
        end
        jinit = jend + 1
    end
    return true
end
