"""
    sign_cadlag(x::N)::N where {N<:Real}

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
julia> import LazySets.sign_cadlag

julia> sign_cadlag.([-0.6, 1.3, 0.0])
3-element Array{Float64,1}:
 -1.0
  1.0
  1.0
```
"""
function sign_cadlag(x::N)::N where {N<:Real}
    return x < zero(x) ? -one(x) : one(x)
end

"""
    ispermutation(u::AbstractVector{T}, v::AbstractVector{T})::Bool where T

Check that two vectors contain the same elements up to reordering.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`true` iff the vectors are identical up to reordering.

### Examples

```jldoctest
julia> LazySets.ispermutation([1, 2, 2], [2, 2, 1])
true

julia> LazySets.ispermutation([1, 2, 2], [1, 1, 2])
false

```
"""
function ispermutation(u::AbstractVector{T}, v::AbstractVector{T})::Bool where T
    if length(u) != length(v)
        return false
    end
    occurrence_map = Dict{T, Int}()
    has_duplicates = false
    for e in u
        if e ∉ v
            return false
        end
        if haskey(occurrence_map, e)
            occurrence_map[e] += 1
            has_duplicates = true
        else
            occurrence_map[e] = 1
        end
    end
    if has_duplicates
        for e in v
            if !haskey(occurrence_map, e) || occurrence_map[e] == 0
                return false
            end
            occurrence_map[e] -= 1
        end
    end
    return true
end

"""
    samedir(u::AbstractVector{N},
            v::AbstractVector{N})::Tuple{Bool, Real} where N<:Real

Check whether two vectors point in the same direction.

### Input

- `u` -- first vector
- `v` -- second vector

### Output

`(true, k)` iff the vectors are identical up to a positive scaling factor `k`,
and `(false, 0)` otherwise.


### Examples

```jldoctest
julia> LazySets.samedir([1, 2, 3], [2, 4, 6])
(true, 0.5)

julia> LazySets.samedir([1, 2, 3], [3, 2, 1])
(false, 0)

julia> LazySets.samedir([1, 2, 3], [-1, -2, -3])
(false, 0)

```
"""
function samedir(u::AbstractVector{N},
                 v::AbstractVector{N}
                )::Tuple{Bool, Real} where N<:Real
    @assert length(u) == length(v) "wrong dimension"
    no_factor = true
    factor = 0
    @inbounds for i in 1:length(u)
        if u[i] == 0
            if v[i] != 0
                return (false, 0)
            end
            continue
        elseif v[i] == 0
            return (false, 0)
        end
        if no_factor
            no_factor = false
            factor = u[i] / v[i]
            if factor < 0
                return (false, 0)
            end
        elseif factor != u[i] / v[i]
            return (false, 0)
        end
    end
    if no_factor
        # both vectors are zero
        return (true, 0)
    end
    return (true, factor)
end

"""
    nonzero_indices(v::AbstractVector{N})::Vector{Int} where N<:Real

Return the indices in which a vector is non-zero.

### Input

- `v` -- vector

### Output

A vector of ascending indices `i` such that the vector is non-zero in dimension
`i`.
"""
function nonzero_indices(v::AbstractVector{N})::Vector{Int} where N<:Real
    n = length(v)
    res = Vector{Int}()
    sizehint!(res, n)
    for i in 1:n
        if v[i] != zero(N)
            push!(res, i)
        end
    end
    return res
end

function nonzero_indices(v::SparseVector{N})::Vector{Int} where N<:Real
    return x.nzind
end

"""
    @neutral(SET, NEUT)

Create functions to make a lazy set operation commutative with a given neutral
element set type.

### Input

- `SET`  -- lazy set operation type
- `NEUT` -- set type for neutral element

### Output

Nothing.

### Notes

This macro generates four functions (possibly two more if `@absorbing` has been
used in advance) (possibly two or four more if `@declare_array_version` has been
used in advance).

### Examples

`@neutral(MinkowskiSum, N)` creates at least the following functions:
* `neutral(::MinkowskiSum) = N`
* `MinkowskiSum(X, N) = X`
* `MinkowskiSum(N, X) = X`
* `MinkowskiSum(N, N) = N`
"""
macro neutral(SET, NEUT)
    @eval begin
        # create function to obtain the neutral element
        function neutral(::Type{$SET})
            return $NEUT
        end

        # create functions to declare the neutral element
        function $SET(X::LazySet{N}, ::$NEUT{N}) where {N<:Real}
            return X
        end
        function $SET(::$NEUT{N}, X::LazySet{N}) where {N<:Real}
            return X
        end
        function $SET(Y::$NEUT{N}, ::$NEUT{N}) where {N<:Real}
            return Y
        end

        # if the absorbing element has already been defined, create combinations
        if isdefined(@__MODULE__, :absorbing) && hasmethod(absorbing, (Type{$SET},))
            @eval(@neutral_absorbing($(esc(SET)),
                                     $(esc(NEUT)),
                                     absorbing($(esc(SET)))))
        end

        # if the array set type has already been defined, create combinations
        if isdefined(@__MODULE__, :array_constructor) &&
                hasmethod(array_constructor, (Type{$SET},))
            @eval(@array_neutral($(esc(SET)),
                                 $(esc(NEUT)),
                                 array_constructor($(esc(SET)))))
        elseif isdefined(@__MODULE__, :is_array_constructor) &&
                hasmethod(is_array_constructor, (Type{$SET},))
            @eval(@array_neutral($(esc(SET)),
                                 $(esc(NEUT)),
                                 $(esc(SET))))
        end
    end
    return nothing
end

"""
    @absorbing(SET, ABS)

Create functions to make a lazy set operation commutative with a given absorbing
element set type.

### Input

- `SET` -- lazy set operation type
- `ABS` -- set type for absorbing element

### Output

Nothing.

### Notes

This macro generates four functions (possibly two more if `@neutral` has been
used in advance) (possibly two or four more if `@declare_array_version` has been
used in advance).

### Examples

`@absorbing(MinkowskiSum, A)` creates at least the following functions:
* `absorbing(::MinkowskiSum) = A`
* `MinkowskiSum(X, A) = A`
* `MinkowskiSum(A, X) = A`
* `MinkowskiSum(A, A) = A`
"""
macro absorbing(SET, ABS)
    @eval begin
        # create function to obtain the absorbing element
        function absorbing(::Type{$SET})
            return $ABS
        end

        # create functions to declare the absorbing element
        function $SET(::LazySet{N}, Y::$ABS{N}) where {N<:Real}
            return Y
        end
        function $SET(Y::$ABS{N}, ::LazySet{N}) where {N<:Real}
            return Y
        end
        function $SET(Y::$ABS{N}, ::$ABS{N}) where {N<:Real}
            return Y
        end

        # if the neutral element has already been defined, create combinations
        if isdefined(@__MODULE__, :neutral) && hasmethod(neutral, (Type{$SET},))
            @eval(@neutral_absorbing($(esc(SET)),
                                     neutral($(esc(SET))),
                                     $(esc(ABS))))
        end

        # if the array set type has already been defined, create combinations
        if isdefined(@__MODULE__, :array_constructor) &&
                hasmethod(array_constructor, (Type{$SET},))
            @eval(@array_absorbing($(esc(SET)),
                                   $(esc(ABS)),
                                   array_constructor($(esc(SET)))))
        elseif isdefined(@__MODULE__, :is_array_constructor) &&
                hasmethod(is_array_constructor, (Type{$SET},))
            @eval(@array_absorbing($(esc(SET)),
                                   $(esc(ABS)),
                                   $(esc(SET))))
        end
    end
    return nothing
end

"""
    @declare_array_version(SET, SETARR)

Create functions to connect a lazy set operation with its array set type.

### Input

- `SET`    -- lazy set operation type
- `SETARR` -- array set type

### Output

Nothing.

### Notes

This macro generates eight functions (and possibly up to eight more if
`@neutral`/`@absorbing` has been used in advance for the base and/or array set
type).

### Examples

`@declare_array_version(MinkowskiSum, MinkowskiSumArray)` creates at least the
following functions:
* `array_constructor(::MinkowskiSum) = MinkowskiSumArray`
* `is_array_constructor(::MinkowskiSumArray) = true`
* `MinkowskiSum!(X, Y)`
* `MinkowskiSum!(X, arr)`
* `MinkowskiSum!(arr, X)`
* `MinkowskiSum!(arr1, arr2)`
"""
macro declare_array_version(SET, SETARR)
    _SET! = Symbol(string(SET), '!')
    @eval begin
        # create function to obtain the array version
        function array_constructor(::Type{$SET})
            return $SETARR
        end

        # create function to check that this is an array version
        function is_array_constructor(::Type{$SETARR})
            return true
        end

        # create in-place modification functions for array version
        function $_SET!(X::LazySet{N}, Y::LazySet{N}) where {N<:Real}
            # no array type: just use the lazy operation
            return $SET(X, Y)
        end
        function $_SET!(X::LazySet{N}, arr::$SETARR{N}) where {N<:Real}
            push!(array(arr), X)
            return arr
        end
        function $_SET!(arr::$SETARR{N}, X::LazySet{N}) where {N<:Real}
            push!(array(arr), X)
            return arr
        end
        function $_SET!(arr1::$SETARR{N}, arr2::$SETARR{N}) where {N<:Real}
            append!(array(arr1), array(arr2))
            return arr1
        end
        # handle method ambiguities with neutral elements
        if isdefined(@__MODULE__, :neutral) && hasmethod(neutral, (Type{$SET},))
            @eval(@array_neutral($(esc(SET)),
                                 neutral($(esc(SET))),
                                 $(esc(SETARR))))
        end
        if isdefined(@__MODULE__, :neutral) && hasmethod(neutral, (Type{$SETARR},))
            @eval(@array_neutral($(esc(SETARR)),
                                 neutral($(esc(SETARR))),
                                 $(esc(SETARR))))
        end
        # handle method ambiguities with absorbing elements
        if isdefined(@__MODULE__, :absorbing) && hasmethod(absorbing, (Type{$SET},))
            @eval(@array_absorbing($(esc(SET)),
                                   absorbing($(esc(SET))),
                                   $(esc(SETARR))))
        end
        if isdefined(@__MODULE__, :absorbing) && hasmethod(absorbing, (Type{$SETARR},))
            @eval(@array_absorbing($(esc(SETARR)),
                                   absorbing($(esc(SETARR))),
                                   $(esc(SETARR))))
        end
    end
    return nothing
end


"""
    @neutral_absorbing(SET, NEUT, ABS)

Create two functions to avoid method ambiguties for a lazy set operation with
respect to neutral and absorbing element set types.

### Input

- `SET`  -- lazy set operation type
- `NEUT` -- set type for neutral element
- `ABS`  -- set type for absorbing element

### Output

A quoted expression containing the function definitions.

### Examples

`@neutral_absorbing(MinkowskiSum, N, A)` creates the following functions as
quoted expressions:
* `MinkowskiSum(N, A) = A`
* `MinkowskiSum(A, N) = A`
"""
macro neutral_absorbing(SET, NEUT, ABS)
    return quote
        function $SET(::$NEUT{N}, Y::$ABS{N}) where {N<:Real}
            return Y
        end
        function $SET(Y::$ABS{N}, ::$NEUT{N}) where {N<:Real}
            return Y
        end
    end
end

"""
    @array_neutral(FUN, NEUT, SETARR)

Create two functions to avoid method ambiguities for a lazy set operation with
respect to the neutral element set type and the array set type.

### Input

- `FUN`     -- function name
- `NEUT`    -- set type for neutral element
- `SETARR`  -- array set type

### Output

A quoted expression containing the function definitions.

### Examples

`@array_neutral(MinkowskiSum, N, ARR)` creates the following functions as
quoted expressions:
* `MinkowskiSum(N, ARR) = ARR`
* `MinkowskiSum(ARR, N) = ARR`
"""
macro array_neutral(FUN, NEUT, SETARR)
    return quote
        function $FUN(::$NEUT{N}, X::$SETARR{N}) where {N<:Real}
            return X
        end
        function $FUN(X::$SETARR{N}, ::$NEUT{N}) where {N<:Real}
            return X
        end
    end
end

"""
    @array_absorbing(FUN, ABS, SETARR)

Create two functions to avoid method ambiguities for a lazy set operation with
respect to the absorbing element set type and the array set type.

### Input

- `FUN`     -- function name
- `ABS`     -- set type for absorbing element
- `SETARR`  -- array set type

### Output

A quoted expression containing the function definitions.

### Examples

`@array_absorbing(MinkowskiSum, ABS, ARR)` creates the following functions as
quoted expressions:
* `MinkowskiSum(ABS, ARR) = ABS`
* `MinkowskiSum(ARR, ABS) = ABS`
"""
macro array_absorbing(FUN, ABS, SETARR)
    return quote
        function $FUN(Y::$ABS{N}, ::$SETARR{N}) where {N<:Real}
            return Y
        end
        function $FUN(::$SETARR{N}, Y::$ABS{N}) where {N<:Real}
            return Y
        end
    end
end
