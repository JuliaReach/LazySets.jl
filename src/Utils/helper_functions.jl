"""
    sign_cadlag(x::N) where {N<:Real}

This function works like the sign function but is ``1`` for input ``0``.

### Input

- `x` -- real scalar

### Output

``1`` if ``x ≥ 0``, ``-1`` otherwise.

### Notes

This is the sign function right-continuous at zero (see
[càdlàg function](https://en.wikipedia.org/wiki/C%C3%A0dl%C3%A0g)).
It can be used with vector-valued arguments via the dot operator.

### Examples

```jldoctest
julia> LazySets.sign_cadlag.([-0.6, 1.3, 0.0])
3-element Array{Float64,1}:
 -1.0
  1.0
  1.0
```
"""
function sign_cadlag(x::N) where {N<:Real}
    return x < zero(x) ? -one(x) : one(x)
end

"""
    substitute(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}

Apply a substitution to a given vector.

### Input

- `substitution` -- substitution (a mapping from an index to a new value)
- `x`            -- vector

### Output

A fresh vector corresponding to `x` after `substitution` was applied.
"""
function substitute(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}
    return substitute!(substitution, copy(x))
end

"""
    substitute!(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}

Apply a substitution to a given vector.

### Input

- `substitution` -- substitution (a mapping from an index to a new value)
- `x`            -- vector (modified in this function)

### Output

The same (but see the Notes below) vector `x` but after `substitution` was
applied.

### Notes

The vector `x` is modified in-place if it has type `Vector` or `SparseVector`.
Otherwise, we first create a new `Vector` from it.
"""
function substitute!(substitution::Dict{Int, T}, x::AbstractVector{T}) where {T}
    return substitute!(Vector(x), substitution)
end

function substitute!(substitution::Dict{Int, T},
                     x::Union{Vector{T}, SparseVector{T}}) where {T}
    for (index, value) in substitution
        x[index] = value
    end
    return x
end

"""
    reseed(rng::AbstractRNG, seed::Union{Int, Nothing})

Reset the RNG seed if the seed argument is a number.

### Input

- `rng`  -- random number generator
- `seed` -- seed for reseeding

### Output

The input RNG if the seed is `nothing`, and a reseeded RNG otherwise.
"""
function reseed(rng::AbstractRNG, seed::Union{Int, Nothing})
    if seed != nothing
        return Random.seed!(rng, seed)
    end
    return rng
end

"""
    StrictlyIncreasingIndices

Iterator over the vectors of `m` strictly increasing indices from 1 to `n`.

### Fields

- `n` -- size of the index domain
- `m` -- number of indices to choose (resp. length of the vectors)

### Notes

The vectors are modified in-place.

The iterator ranges over ``\\binom{n}{m}`` (`n` choose `m`) possible vectors.

This implementation results in a lexicographic order with the last index growing
first.

### Examples

```jldoctest
julia> for v in LazySets.StrictlyIncreasingIndices(4, 2)
           println(v)
       end
[1, 2]
[1, 3]
[1, 4]
[2, 3]
[2, 4]
[3, 4]
```
"""
struct StrictlyIncreasingIndices
    n::Int
    m::Int

    function StrictlyIncreasingIndices(n::Int, m::Int)
        @assert n >= m > 0 "require n >= m > 0"
        new(n, m)
    end
end

Base.eltype(::Type{StrictlyIncreasingIndices}) = Vector{Int}
Base.length(sii::StrictlyIncreasingIndices) = binomial(sii.n, sii.m)

# initialization
function Base.iterate(sii::StrictlyIncreasingIndices)
    v = [1:sii.m;]
    return (v, v)
end

# normal iteration
function Base.iterate(sii::StrictlyIncreasingIndices, state::AbstractVector{Int})
    v = state
    i = sii.m
    diff = sii.n
    if i == diff
        return nothing
    end
    while v[i] == diff
        i -= 1
        diff -= 1
    end
    # update vector
    v[i] += 1
    for j in i+1:sii.m
        v[j] = v[j-1] + 1
    end
    # detect termination: first index has maximum value
    if i == 1 && v[1] == (sii.n - sii.m + 1)
        return (v, nothing)
    end
    return (v, v)
end

# termination
function Base.iterate(sii::StrictlyIncreasingIndices, state::Nothing)
    return nothing
end

"""
    subtypes(interface, concrete::Bool)

Return the concrete subtypes of a given interface.

### Input

- `interface` -- an abstract type, usually a set interface
- `concrete`  -- if `true`, seek further the inner abstract subtypes of the given
                 interface, otherwise return only the direct subtypes of `interface`

### Output

A list with the subtypes of the abstract type `interface`, sorted alphabetically.

### Examples

Consider the `AbstractPolytope` interface. If we include the abstract subtypes
of this interface,

```jldoctest subtypes
julia> using LazySets: subtypes

julia> subtypes(AbstractPolytope, false)
4-element Array{Any,1}:
 AbstractCentrallySymmetricPolytope
 AbstractPolygon
 HPolytope
 VPolytope
```

We can use this function to obtain the concrete subtypes of
`AbstractCentrallySymmetricPolytope` and `AbstractPolygon` (further until all
concrete types are obtained), using the `concrete` flag:

```jldoctest subtypes
julia> subtypes(AbstractPolytope, true)
15-element Array{Type,1}:
 Ball1
 BallInf
 HParallelotope
 HPolygon
 HPolygonOpt
 HPolytope
 Hyperrectangle
 Interval
 LineSegment
 Singleton
 SymmetricIntervalHull
 VPolygon
 VPolytope
 ZeroSet
 Zonotope
```
"""
function subtypes(interface, concrete::Bool)

    subtypes_to_test = subtypes(interface)

    # do not seek the concrete subtypes further
    if !concrete
        return sort(subtypes_to_test, by=string)
    end

    result = Vector{Type}()
    i = 0
    while i < length(subtypes_to_test)
        i += 1
        subtype = subtypes_to_test[i]
        new_subtypes = subtypes(subtype)
        if isempty(new_subtypes)
            # base type found
            push!(result, subtype)
        else
            # yet another interface layer
            append!(subtypes_to_test, new_subtypes)
        end
    end
    return sort(result, by=string)
end

"""
    implementing_sets(op::Function;
                      signature::Tuple{Vector{Type}, Int}=(Type[], 1),
                      type_args=Float64, binary::Bool=false)

Compute a dictionary containing information about availability of (unary or
binary) concrete set operations.

### Input

- `op`        -- set operation (respectively its `Function` object)
- `signature` -- (optional, default: `Type[]`) the type signature of the
                 function without the `LazySet` type(s) (see also the `index`
                 option and the `Examples` section below)
- `index`     -- (optional, default: `1`) index of the set type in the signature
                 in the unary case (see the `binary` option)
- `type_args` -- (optional, default: `Float64`) type arguments added to the
                 `LazySet`(s) when searching for available methods; valid
                 inputs are a type or `nothing`, and in the unary case (see the
                 `binary` option) it can also be a list of types
- `binary`    -- (optional, default: `false`) flag indicating whether `op` is a
                 binary function (`true`) or a unary function (`false`)

### Output

A dictionary with three keys each mapping to a list:
- `"available"` -- This list contains all set types such that there exists an
                   implementation of `op`.
- `"missing"`   -- This list contains all set types such that there does not
                   exist an implementation of `op`. Note that this is the
                   complement of the `"available"` list.
- `"specific"`  -- This list contains all set types such that there exists a
                   type-specific implementation. Note that those set types also
                   occur in the `"available"` list.

In the unary case, the lists contain set types.
In the binary case, the lists contain pairs of set types.

### Examples

```jldoctests
julia> using LazySets: implementing_sets

julia> dict = implementing_sets(tovrep);

julia> dict["available"]  # tovrep is only available for polyhedral set types
6-element Array{Type,1}:
 HPolygon
 HPolygonOpt
 HPolyhedron
 HPolytope
 VPolygon
 VPolytope

julia> dict = implementing_sets(σ; signature=Type[AbstractVector{Float64}], index=2);

julia> dict["missing"]  # every set type implements function σ
0-element Array{Type,1}

julia> N = Rational{Int};  # restriction of the number type

julia> dict = implementing_sets(σ; signature=Type[AbstractVector{N}], index=2, type_args=N);

julia> dict["missing"]  # some set types are not available with number type N
4-element Array{Type,1}:
 Ball2
 Ballp
 Bloating
 Ellipsoid

julia> dict = LazySets.implementing_sets(convex_hull; binary=true);  # binary case

julia> (HPolytope, HPolytope) ∈ dict["available"]  # dict contains pairs now
true
```
"""
function implementing_sets(op::Function;
                           signature::Vector{Type}=Type[],
                           index::Int=1,
                           type_args=Float64,
                           binary::Bool=false)
    function get_unary_dict()
        return Dict{String, Vector{Type}}("available" => Vector{Type}(),
                                          "missing" => Vector{Type}(),
                                          "specific" => Vector{Type}())
    end

    if !binary
        # unary case
        dict = get_unary_dict()
        _implementing_sets_unary!(dict, op, signature, index, type_args)
        return dict
    end

    # binary case by considering all set types as second argument
    binary_dict = Dict{String, Vector{Tuple{Type, Type}}}(
        "available" => Vector{Tuple{Type, Type}}(),
        "missing" => Vector{Tuple{Type, Type}}(),
        "specific" => Vector{Tuple{Type, Type}}())
    signature_copy = copy(signature)
    insert!(signature_copy, 1, LazySet)  # insert a dummy type
    for set_type in subtypes(LazySet, true)
        augmented_set_type = set_type
        try
            if type_args isa Type
                augmented_set_type = set_type{type_args}
            else
                @assert type_args == nothing "for binary functions, " *
                    "`type_args` must not be a list"
            end
        catch e
            # type combination not available
            push!(binary_dict["missing"], (set_type, LazySet))
            continue
        end
        signature_copy[1] = augmented_set_type  # replace dummy type
        unary_dict = get_unary_dict()
        _implementing_sets_unary!(unary_dict, op, signature_copy, 1, type_args)

        # create pairs
        for key in ["available", "missing", "specific"]
            i = 1
            for other_set_type in unary_dict[key]
                push!(binary_dict[key], (set_type, other_set_type))
                i += 1
            end
        end
    end
    return binary_dict
end

function _implementing_sets_unary!(dict, op, signature, index, type_args)
    function create_type_tuple(given_type)
        if isempty(signature)
            return Tuple{given_type}
        end
        # create dynamic tuple
        new_signature = copy(signature)
        insert!(new_signature, index, given_type)
        return Tuple{new_signature...}
    end

    for set_type in subtypes(LazySet, true)
        augmented_set_type = set_type
        try
            if type_args isa Type
                augmented_set_type = set_type{type_args}
            elseif type_args != nothing
                augmented_set_type = set_type{type_args...}
            end
        catch e
            # type combination not available
            push!(dict["missing"], set_type)
            continue
        end
        tuple_set_type = create_type_tuple(augmented_set_type)
        if !hasmethod(op, tuple_set_type)
            # implementation not available
            push!(dict["missing"], set_type)
            continue
        end
        # implementation available
        push!(dict["available"], set_type)

        # check if there is a type-specific method
        super_type = supertype(augmented_set_type)
        tuple_super_type = create_type_tuple(super_type)
        if !hasmethod(op, tuple_super_type)
            continue
        end
        method_set_type = which(op, tuple_set_type)
        method_super_type = which(op, tuple_super_type)
        if method_set_type != method_super_type
            push!(dict["specific"], set_type)
        end
    end
end

# check that the given coordinate i can be used to index an arbitrary element in the set X
@inline function _check_bounds(X, i)
    1 <= i <= dim(X) || throw(ArgumentError("there is no index at coordinate $i, since the set is of dimension $(dim(X))"))
end
