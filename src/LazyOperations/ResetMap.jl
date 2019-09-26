import Base: isempty

export ResetMap,
       get_A,
       get_b

"""
    ResetMap{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents a lazy reset map.
A reset map is a special case of an affine map ``A x + b, x ∈ X`` where the
linear map ``A`` is the identity matrix with zero entries in all reset
dimensions, and the translation vector ``b`` is zero in all other dimensions.

### Fields

- `X`      -- convex set
- `resets` -- resets (a mapping from an index to a new value)

### Examples

```jldoctest resetmap
julia> X = BallInf([2.0, 2.0, 2.0], 1.0);

julia> r = Dict(1 => 4.0, 3 => 0.0);

julia> rm = ResetMap(X, r);

```

Here `rm` modifies the set `X` such that `x1` is reset to 4 and `x3` is reset to
0, while `x2` is not modified.
Hence `rm` is equivalent to the set
`Hyperrectangle([4.0, 2.0, 0.0], [0.0, 1.0, 0.0])`, i.e., an axis-aligned line
segment embedded in 3D.

The corresponding affine map ``A x + b`` would be:

```math
    \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix} x +
    \begin{pmatrix} 4 & 0 & 0 \end{pmatrix}
```

Use the function `get_A` (resp. `get_b`) to create the matrix `A` (resp. vector
`b`) corresponding to a given reset map.

```jldoctest resetmap
julia> get_A(rm)
3×3 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:
 0.0   ⋅    ⋅
  ⋅   1.0   ⋅
  ⋅    ⋅   0.0

julia> get_b(rm)
3-element SparseArrays.SparseVector{Float64,Int64} with 1 stored entry:
  [1]  =  4.0
```

The application of a `ResetMap` to a `ZeroSet` or an `EmptySet` is simplified
automatically.

```jldoctest resetmap
julia> ResetMap(ZeroSet(3), r)
Singleton{Float64,SparseArrays.SparseVector{Float64,Int64}}(  [1]  =  4.0)

julia> ResetMap(EmptySet(), r)
EmptySet{Float64}()
```

The (in this case unique) support vector of `rm` in direction `ones(3)` is:

```jldoctest resetmap
julia> σ(ones(3), rm)
3-element Array{Float64,1}:
 4.0
 3.0
 0.0
```
"""
struct ResetMap{N<:Real, S<:LazySet{N}} <: LazySet{N}
    X::S
    resets::Dict{Int, N}
end

isoperation(::Type{ResetMap}) = true

# ZeroSet is "almost absorbing" for the linear map (only the dimension changes)
# such that only the translation vector remains
function ResetMap(Z::ZeroSet{N}, resets::Dict{Int, N}) where {N<:Real}
    return Singleton(_get_b_from_dictionary(resets, dim(Z)))
end

# EmptySet is absorbing for ResetMap
function ResetMap(∅::EmptySet{N}, resets::Dict{Int, N}) where {N<:Real}
    return ∅
end

"""
    get_A(rm::ResetMap{N}) where {N<:Real}

Return the ``A`` matrix of the affine map ``A x + b, x ∈ X`` represented by a
reset map.

### Input

- `rm` -- reset map

### Output

The (diagonal) matrix for the affine map ``A x + b, x ∈ X`` represented by the
reset map.

### Algorithm

We construct the identity matrix and set all entries in the reset dimensions to
zero.
"""
function get_A(rm::ResetMap{N}) where {N<:Real}
    n = dim(rm)
    v = ones(N, n)
    for i in keys(rm.resets)
        v[i] = zero(N)
    end
    return Diagonal(v)
end

"""
    get_b(rm::ResetMap{N}) where {N<:Real}

Return the ``b`` vector of the affine map ``A x + b, x ∈ X`` represented by a
reset map.

### Input

- `rm` -- reset map

### Output

The (sparse) vector for the affine map ``A x + b, x ∈ X`` represented by the
reset map.
The vector contains the reset value for all reset dimensions, and is zero for
all other dimensions.
"""
function get_b(rm::ResetMap{N}) where {N<:Real}
    return _get_b_from_dictionary(rm.resets, dim(rm))
end

function _get_b_from_dictionary(dict::Dict{Int, N}, n::Int) where {N<:Real}
    b = sparsevec(Int[], N[], n)
    for (i, val) in dict
        b[i] = val
    end
    return b
end


# --- LazySet interface functions ---


"""
    dim(rm::ResetMap)

Return the dimension of a reset map.

### Input

- `rm` -- reset map

### Output

The dimension of a reset map.
"""
function dim(rm::ResetMap)::Int
    return dim(rm.X)
end

"""
    σ(d::AbstractVector{N}, rm::ResetMap{N}) where {N<:Real}

Return the support vector of a reset map.

### Input

- `d`  -- direction
- `rm` -- reset map

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
function σ(d::AbstractVector{N}, rm::ResetMap{N}) where {N<:Real}
    d_reset = copy(d)
    for var in keys(rm.resets)
        d_reset[var] = zero(N)
    end
    return substitute(rm.resets, σ(d_reset, rm.X))
end

"""
    ρ(d::AbstractVector{N}, rm::ResetMap{N}) where {N<:Real}

Return the support function of a reset map.

### Input

- `d`  -- direction
- `rm` -- reset map

### Output

The support function in the given direction.

### Notes

We use the usual dot-product definition, but for unbounded sets we redefine the
product between ``0`` and ``±∞`` as ``0``; Julia returns `NaN` here.

```jldoctest
julia> Inf * 0.0
NaN
```

See the discussion
[here](https://math.stackexchange.com/questions/28940/why-is-infty-cdot-0-not-clearly-equal-to-0).
"""
function ρ(d::AbstractVector{N}, rm::ResetMap{N}) where {N<:Real}
    return dot_zero(d, σ(d, rm))
end

"""
    an_element(rm::ResetMap)

Return some element of a reset map.

### Input

- `rm` -- reset map

### Output

An element in the reset map.
It relies on the `an_element` function of the wrapped set.
"""
function an_element(rm::ResetMap)
    return substitute(rm.resets, an_element(rm.X))
end

"""
    isempty(rm::ResetMap)::Bool

Return if a reset map is empty or not.

### Input

- `rm` -- reset map

### Output

`true` iff the wrapped set is empty.
"""
function isempty(rm::ResetMap)::Bool
    return isempty(rm.X)
end

"""
    constraints_list(rm::ResetMap{N}) where {N<:Real}

Return the list of constraints of a polytopic reset map.

### Input

- `rm` -- reset map of a polytope

### Output

The list of constraints of the reset map.

### Notes

We assume that the underlying set `X` is a polytope, i.e., is bounded and offers
a method `constraints_list(X)`.

### Algorithm

We fall back to `constraints_list` of a `LinearMap` of the `A`-matrix in the
affine-map view of a reset map.
Each reset dimension ``i`` is projected to zero, expressed by two constraints
for each reset dimension.
Then it remains to shift these constraints to the new value.

For instance, if the dimension ``5`` was reset to ``4``, then there will be
constraints ``x₅ ≤ 0`` and ``-x₅ ≤ 0``.
We then modify the right-hand side of these constraints to ``x₅ ≤ 4`` and
``-x₅ ≤ -4``, respectively.
"""
function constraints_list(rm::ResetMap{N}) where {N<:Real}
    # if `vector` has exactly one non-zero entry, return its index
    # otherwise return 0
    function find_unique_nonzero_entry(vector::AbstractVector{N})
        res = 0
        for (i, v) in enumerate(vector)
            if v != zero(N)
                if res != 0
                    # at least two non-zero entries
                    return 0
                else
                    # first non-zero entry so far
                    res = i
                end
            end
        end
        return res
    end

    constraints = copy(constraints_list(LinearMap(get_A(rm), rm.X)))
    for (i, c) in enumerate(constraints)
        constrained_dim = find_unique_nonzero_entry(c.a)
        if constrained_dim > 0  # constraint in only one dimension
            if !haskey(rm.resets, constrained_dim)
                continue  # not a dimension we are interested in
            end
            new_value = rm.resets[constrained_dim]
            if new_value == zero(N)
                @assert c.b == zero(N)
                continue  # a reset to 0 needs not create a new constraint
            end
            if c.a[constrained_dim] < zero(N)
                # change sign for lower bound
                new_value = -new_value
            end
            constraints[i] = HalfSpace(c.a, new_value)
        end
    end
    return constraints
end

"""
    constraints_list(rm::ResetMap{N, S}) where
        {N<:Real, S<:AbstractHyperrectangle}

Return the list of constraints of a hyperrectangular reset map.

### Input

- `rm` -- reset map of a hyperrectangular set

### Output

The list of constraints of the reset map.

### Algorithm

We iterate through all dimensions.
If there is a reset, we construct the corresponding (flat) constraints.
Otherwise, we construct the corresponding constraints of the underlying set.
"""
function constraints_list(rm::ResetMap{N, S}
                         ) where {N<:Real, S<:AbstractHyperrectangle}
    H = rm.X
    n = dim(H)
    constraints = Vector{LinearConstraint{N}}(undef, 2*n)
    j = 1
    for i in 1:n
        ei = SingleEntryVector(i, n, one(N))
        if haskey(rm.resets, i)
            # reset dimension => add flat constraints
            v = rm.resets[i]
            constraints[j] = HalfSpace(ei, v)
            constraints[j+1] = HalfSpace(-ei, -v)
        else
            # non-reset dimension => use the hyperrectangle's constraints
            constraints[j] = HalfSpace(ei, high(H, i))
            constraints[j+1] = HalfSpace(-ei, -low(H, i))
        end
        j += 2
    end
    return constraints
end
