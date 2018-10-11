import Base:*, ∈

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

@static if VERSION < v"0.7-"
    # convenience constructor without type parameter
    CartesianProduct(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} =
        CartesianProduct{N, S1, S2}(X, Y)
end

# constructor from an array
CartesianProduct(Xarr::Vector{S}) where {N<:Real, S<:LazySet{N}} =
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
    σ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}

Return the support vector of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function σ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}
    n1 = dim(cp.X)
    return [σ(d[1:n1], cp.X); σ(d[n1+1:length(d)], cp.Y)]
end

"""
    ρ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}

Return the support function of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function ρ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}
    n1 = dim(cp.X)
    return ρ(d[1:n1], cp.X) + ρ(d[n1+1:length(d)], cp.Y)
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

    n1 = dim(cp.X)
    return ∈(view(x, 1:n1), cp.X) &&
           ∈(view(x, n1+1:length(x)), cp.Y)
end

"""
    constraints_list(cp::CartesianProduct{N})::Vector{LinearConstraint{N}} where N<:Real

Return the list of constraints of a (polytopic) Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

A list of constraints.
"""
function constraints_list(cp::CartesianProduct{N})::Vector{LinearConstraint{N}} where N<:Real
    return constraints_list(CartesianProductArray([cp.X, cp.Y]))
end

"""
    vertices_list(cp::CartesianProduct{N})::Vector{Vector{N}} where N<:Real

Return the list of vertices of a (polytopic) Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
function vertices_list(cp::CartesianProduct{N})::Vector{Vector{N}} where N<:Real
    # collect low-dimensional vertices lists
    vlist_low = (vertices_list(cp.X), vertices_list(cp.Y))

    # create high-dimensional vertices list
    vlist = Vector{Vector{N}}()
    m = length(vlist_low[1]) * length(vlist_low[2])
    sizehint!(vlist, m)
    for v1 in vlist_low[1]
        for v2 in vlist_low[2]
            push!(vlist, vcat(v1, v2))
        end
    end

    return vlist
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

@static if VERSION < v"0.7-"
    # convenience constructor without type parameter
    CartesianProductArray(arr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
        CartesianProductArray{N, S}(arr)
end

# constructor for an empty product with optional size hint and numeric type
function CartesianProductArray(n::Int=0, N::Type=Float64)::CartesianProductArray
    arr = Vector{LazySet{N}}()
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
    return length(cpa.array) == 0 ? 0 : sum([dim(Xi) for Xi in cpa.array])
end

"""
    σ(d::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}

Support vector of a Cartesian product array.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the product sets.
"""
function σ(d::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}
    svec = similar(d)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        svec[i0:i1] = σ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return svec
end

"""
    ρ(d::AbstractVector{N}, cp::CartesianProductArray{N}) where {N<:Real}

Return the support function of a Cartesian product array.

### Input

- `d`  -- direction
- `cpa` -- Cartesian product array

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function ρ(d::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}
    sfun = zero(N)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        sfun += ρ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return sfun
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

    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if !∈(x[i0:i1], Xi)
            return false
        end
        i0 = i1 + 1
    end
    return true
end

"""
    constraints_list(cpa::CartesianProductArray{N})::Vector{LinearConstraint{N}} where N<:Real

Return the list of constraints of a (polytopic) Cartesian product.

### Input

- `cpa` -- Cartesian product

### Output

A list of constraints.

"""
function constraints_list(cpa::CartesianProductArray{N})::Vector{LinearConstraint{N}} where N<:Real
    # collect low-dimensional constraints lists
    c_array = array(cpa)
    clist = Vector{LinearConstraint{N}}()
    sizehint!(clist, sum(dim(s_low) for s_low in c_array))
    prev_step = 1
    # create high-dimensional constraints list
    for c_low in c_array
        if c_low isa LinearConstraint
            indices = prev_step : (dim(c_low) + prev_step - 1)
            new_constr = LinearConstraint(sparsevec(indices, c_low.a), c_low.b)
            push!(clist, new_constr)
            prev_step += dim(c_low)
        else
            c_low_list = constraints_list(c_low)
            if !isempty(c_low_list)
                indices = prev_step : (dim(c_low_list[1]) + prev_step - 1)
            end
            for constr in c_low_list
                new_constr = LinearConstraint(sparsevec(indices, constr.a), constr.b)
                push!(clist, new_constr)
            end
            prev_step += dim(c_low_list[1])
        end
    end

    return clist
end

"""
    vertices_list(cpa::CartesianProductArray{N})::Vector{Vector{N}} where N<:Real

Return the list of vertices of a (polytopic) Cartesian product.

### Input

- `cpa` -- Cartesian product

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
function vertices_list(cpa::CartesianProductArray{N})::Vector{Vector{N}} where N<:Real
    # collect low-dimensional vertices lists
    vlist_low = [vertices_list(X) for X in array(cpa)]

    # create high-dimensional vertices list
    indices_max = [length(vl) for vl in vlist_low]
    m = prod(indices_max)
    vlist = Vector{Vector{N}}(undef, m)
    indices = ones(Int, length(vlist_low))
    v = zeros(N, dim(cpa))
    dim_start_j = 1
    for vl in vlist_low
        v_low = vl[1]
        v[dim_start_j:dim_start_j+length(v_low)-1] = v_low
        dim_start_j += length(v_low)
    end
    i = 1
    j = 1
    # iterate through all index combinations
    while true
        indices[1] = 0
        while indices[1] < indices_max[1]
            indices[1] += 1
            v_low = vlist_low[1][indices[1]]
            v[1:length(v_low)] = v_low
            vlist[i] = copy(v)
            i += 1
        end
        if i > m
            break
        end
        j = 1
        dim_start_j = 1
        while indices[j] == indices_max[j]
            indices[j] = 1
            v_low = vlist_low[j][1]
            v[dim_start_j:dim_start_j+length(v_low)-1] = v_low
            dim_start_j += length(v_low)
            j += 1
        end
        indices[j] += 1
        v_low = vlist_low[j][indices[j]]
        v[dim_start_j:dim_start_j+length(v_low)-1] = v_low
    end

    return vlist
end
