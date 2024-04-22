using LazySets.Arrays: find_unique_nonzero_entry

export ResetMap

"""
    ResetMap{N, S<:LazySet{N}} <: AbstractAffineMap{N, S}

Type that represents a lazy reset map.
A reset map is a special case of an affine map ``A x + b, x ∈ X`` where the
linear map ``A`` is the identity matrix with zero entries in all reset
dimensions, and the translation vector ``b`` is zero in all other dimensions.

### Fields

- `X`      -- set
- `resets` -- resets (a mapping from an index to a new value)

### Notes

The reset map preserves convexity: if `X` is convex, then any reset map of `X`
is convex as well.

### Examples

```jldoctest resetmap
julia> X = BallInf([2.0, 2.0, 2.0], 1.0);

julia> r = Dict(1 => 4.0, 3 => 0.0);

julia> rm = ResetMap(X, r);

```

Here `rm` modifies the set `X` such that `x1` is reset to 4 and `x3` is reset to
0, while `x2` is not modified.
Hence `rm` is equivalent to the set
`VPolytope([[4.0, 1.0, 0.0], [4.0, 3.0, 0.0]])`, i.e., an axis-aligned line
segment embedded in 3D.

The corresponding affine map ``A x + b`` would be:

```math
    \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix} x +
    \begin{pmatrix} 4 & 0 & 0 \end{pmatrix}
```

Use the function `matrix` (resp. `vector`) to create the matrix `A` (resp.
vector `b`) corresponding to a given reset map.

```jldoctest resetmap
julia> matrix(rm)
3×3 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.0   ⋅    ⋅
  ⋅   1.0   ⋅
  ⋅    ⋅   0.0

julia> vector(rm)
3-element SparseArrays.SparseVector{Float64, Int64} with 1 stored entry:
  [1]  =  4.0
```

The application of a `ResetMap` to a `ZeroSet` or an `EmptySet` is simplified
automatically.

```jldoctest resetmap
julia> ResetMap(ZeroSet(3), r)
Singleton{Float64, SparseArrays.SparseVector{Float64, Int64}}(  [1]  =  4.0)

julia> ResetMap(EmptySet(3), r)
∅(3)
```

The (in this case unique) support vector of `rm` in direction `[1, 1, 1]` is:

```jldoctest resetmap
julia> σ(ones(3), rm)
3-element Vector{Float64}:
 4.0
 3.0
 0.0
```
"""
struct ResetMap{N,S<:LazySet{N}} <: AbstractAffineMap{N,S}
    X::S
    resets::Dict{Int,N}
end

isoperationtype(::Type{<:ResetMap}) = true

isconvextype(::Type{ResetMap{N,S}}) where {N,S} = isconvextype(S)

# ZeroSet is "almost absorbing" for the reset map because only the translation
# vector remains
function ResetMap(Z::ZeroSet{N}, resets::Dict{Int,N}) where {N}
    return Singleton(_vector_from_dictionary(resets, dim(Z)))
end

# EmptySet is absorbing for ResetMap
function ResetMap(∅::EmptySet{N}, resets::Dict{Int,N}) where {N}
    return ∅
end

"""
    matrix(rm::ResetMap{N}) where {N}

Return the ``A`` matrix of the affine map ``A x + b, x ∈ X`` represented by a
reset map.

### Input

- `rm` -- reset map

### Output

The (`Diagonal`) matrix for the affine map ``A x + b, x ∈ X`` represented by the
reset map.

### Algorithm

We construct the identity matrix and set all entries in the reset dimensions to
zero.
"""
function matrix(rm::ResetMap{N}) where {N}
    n = dim(rm)
    v = ones(N, n)
    for i in keys(rm.resets)
        v[i] = zero(N)
    end
    return Diagonal(v)
end

"""
    vector(rm::ResetMap)

Return the ``b`` vector of the affine map ``A x + b, x ∈ X`` represented by a
reset map.

### Input

- `rm` -- reset map

### Output

The (sparse) vector for the affine map ``A x + b, x ∈ X`` represented by the
reset map.
The vector contains the reset value for all reset dimensions and is zero for all
other dimensions.
"""
function vector(rm::ResetMap)
    return _vector_from_dictionary(rm.resets, dim(rm))
end

function _vector_from_dictionary(dict::Dict{Int,N}, n::Int) where {N}
    b = sparsevec(Int[], N[], n)
    for (i, val) in dict
        b[i] = val
    end
    return b
end

"""
    set(rm::ResetMap)

Return the set wrapped by a reset map.

### Input

- `rm` -- reset map

### Output

The wrapped set.
"""
function set(rm::ResetMap)
    return rm.X
end

"""
    dim(rm::ResetMap)

Return the dimension of a reset map.

### Input

- `rm` -- reset map

### Output

The ambient dimension of a reset map.
"""
function dim(rm::ResetMap)
    return dim(rm.X)
end

"""
    σ(d::AbstractVector, rm::ResetMap)

Return a support vector of a reset map.

### Input

- `d`  -- direction
- `rm` -- reset map

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
function σ(d::AbstractVector, rm::ResetMap)
    N = promote_type(eltype(d), eltype(rm))
    d_reset = copy(d)
    for var in keys(rm.resets)
        d_reset[var] = zero(N)
    end
    return substitute(rm.resets, σ(d_reset, rm.X))
end

"""
    ρ(d::AbstractVector, rm::ResetMap)

Evaluate the support function of a reset map.

### Input

- `d`  -- direction
- `rm` -- reset map

### Output

The evaluation of the support function in the given direction.

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
function ρ(d::AbstractVector, rm::ResetMap)
    return dot_zero(d, σ(d, rm))
end

"""
    an_element(rm::ResetMap)

Return some element of a reset map.

### Input

- `rm` -- reset map

### Output

An element in the reset map.

### Algorithm

This method relies on the `an_element` implementation for the wrapped set.
"""
function an_element(rm::ResetMap)
    return substitute(rm.resets, an_element(rm.X))
end

function isboundedtype(::Type{<:ResetMap{N,S}}) where {N,S}
    return isboundedtype(S)
end

"""
    constraints_list(rm::ResetMap)

Return a list of constraints of a polyhedral reset map.

### Input

- `rm` -- reset map of a polyhedron

### Output

A list of constraints of the reset map.

### Notes

We assume that the underlying set `rm.X` is a polyhedron, i.e., offers a method
`constraints_list(X)`.

### Algorithm

If the set `rm.X` is hyperrectangular, we iterate through all dimensions.
For each reset we construct the corresponding (flat) constraints, and in the
other dimensions we construct the corresponding constraints of the underlying
set.

For more general sets, we fall back to `constraints_list` of a `LinearMap` of
the `A`-matrix in the affine-map view of a reset map.
Each reset dimension ``i`` is projected to zero, expressed by two constraints
for each reset dimension.
Then it remains to shift these constraints to the new value.

For instance, if the dimension ``5`` was reset to ``4``, then there will be
constraints ``x₅ ≤ 0`` and ``-x₅ ≤ 0``.
We then modify the right-hand side of these constraints to ``x₅ ≤ 4`` and
``-x₅ ≤ -4``, respectively.
"""
function constraints_list(rm::ResetMap)
    constraints = constraints_list(LinearMap(matrix(rm), set(rm)))
    N = eltype(rm)
    for (i, c) in enumerate(constraints)
        constrained_dim = find_unique_nonzero_entry(c.a)
        if constrained_dim > 0  # constrained in only one dimension
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

function constraints_list(rm::ResetMap{N,S}) where {N,S<:AbstractHyperrectangle}
    H = rm.X
    n = dim(H)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2 * n)
    j = 1
    for i in 1:n
        ei = SingleEntryVector(i, n, one(N))
        if haskey(rm.resets, i)
            # reset dimension => add flat constraints
            v = rm.resets[i]
            constraints[j] = HalfSpace(ei, v)
            constraints[j + 1] = HalfSpace(-ei, -v)
        else
            # non-reset dimension => use the hyperrectangle's constraints
            constraints[j] = HalfSpace(ei, high(H, i))
            constraints[j + 1] = HalfSpace(-ei, -low(H, i))
        end
        j += 2
    end
    return constraints
end
