"""
    @neutral(SET, NEUT)

Create methods to make a lazy set operation commutative with a given
neutral-element set type.

### Input

- `SET`  -- set type of lazy operation
- `NEUT` -- set type of neutral element

### Output

Nothing.

### Notes

This macro generates four functions (possibly two more if `@absorbing` has been
used in advance, and possibly two or four more if `@declare_array_version` has
been used in advance).

### Examples

`@neutral(MinkowskiSum, N)` creates at least the following methods:
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
        function $SET(X::LazySet{N}, ::$NEUT{N}) where {N}
            return X
        end
        function $SET(::$NEUT{N}, X::LazySet{N}) where {N}
            return X
        end
        function $SET(Y::$NEUT{N}, ::$NEUT{N}) where {N}
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

Create methods to make a lazy set operation commutative with a given
absorbing-element set type.

### Input

- `SET` -- set type of lazy operation
- `ABS` -- set type of absorbing element

### Output

Nothing.

### Notes

This macro generates four functions (possibly two more if `@neutral` has been
used in advance, and possibly two or four more if `@declare_array_version` has
been used in advance).

### Examples

`@absorbing(MinkowskiSum, A)` creates at least the following methods:
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
        function $SET(::LazySet{N}, Y::$ABS{N}) where {N}
            return Y
        end
        function $SET(Y::$ABS{N}, ::LazySet{N}) where {N}
            return Y
        end
        function $SET(Y::$ABS{N}, ::$ABS{N}) where {N}
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
    @declare_binary_operation(SET)

Create common methods for binary set operations.

### Input

- `SET` -- set type of the lazy operation

### Output

Nothing.

### Notes

This macro generates seven methods. See the example below.

### Examples

`@declare_binary_operation(MinkowskiSum)` creates the following methods:
* `iterate(::MinkowskiSum)`
* `length(::MinkowskiSum)`
* `getindex(::MinkowskiSum, ::Int)`
* `getindex(::MinkowskiSum, ::AbstractVector{Int})`
* `lastindex(::MinkowskiSum)`
* `array(::MinkowskiSum)`
* `is_array_constructor(::Type{MinkowskiSum})`
"""
macro declare_binary_operation(SET)
    @eval begin
        function Base.iterate(X::$SET, state=1)
            if state == 1
                return (first(X), 2)
            elseif state == 2
                return (second(X), 3)
            else
                return nothing
            end
        end

        function Base.length(::$SET)
            return 2
        end

        function Base.lastindex(X::$SET)
            return 2
        end

        function Base.getindex(X::$SET, i::Int)
            if i == 1
                return first(X)
            elseif i == 2
                return second(X)
            end
            throw(ArgumentError("invalid index $i for binary set operation"))
        end

        function Base.getindex(X::$SET, indices::AbstractVector{Int})
            return [X[i] for i in indices]
        end

        function array(X::$SET)
            return [first(X), second(X)]
        end

        function is_array_constructor(::Type{$SET})
            return false
        end

        if hasmethod(concrete_function, (Type{<:$SET},))
            function concretize(X::$SET)
                f = concrete_function($SET)
                return f(concretize(first(X)), concretize(second(X)))
            end
        end
    end
    return nothing
end

"""
    @declare_array_version(SET, SETARR)

Create methods to connect a lazy set operation with its array set type.

### Input

- `SET`    -- set type of lazy operation
- `SETARR` -- set type of array version

### Output

Nothing.

### Notes

This macro generates six methods (and possibly up to eight more if
`@neutral`/`@absorbing` has been used in advance for the base and/or array set
type). See the example below.

### Examples

`@declare_array_version(MinkowskiSum, MinkowskiSumArray)` creates at least the
following methods:
* `array_constructor(::Type{MinkowskiSum}) = MinkowskiSumArray`
* `binary_constructor(::Type{MinkowskiSumArray}) = MinkowskiSum`
* `is_array_constructor(::Type{MinkowskiSumArray}) = true`
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

        # create function to obtain the binary version
        function binary_constructor(::Type{$SETARR})
            return $SET
        end

        # create function to check whether an operation is the array version
        function is_array_constructor(::Type{$SETARR})
            return true
        end

        # create function to flatten a lazy set operation
        function flatten(X::Union{<:$SET,<:$SETARR})
            arr = flatten!([], X, $SET)
            @inbounds if length(arr) == 1
                return arr[1]
            elseif length(arr) == 2
                return $SET(2)
            else
                return $SETARR([Xi for Xi in arr])
            end
        end

        # create in-place modification functions for array version
        function $_SET!(X::LazySet{N}, Y::LazySet{N}) where {N}
            # no array type: just use the lazy operation
            return $SET(X, Y)
        end
        function $_SET!(X::LazySet{N}, arr::$SETARR{N}) where {N}
            push!(array(arr), X)
            return arr
        end
        function $_SET!(arr::$SETARR{N}, X::LazySet{N}) where {N}
            push!(array(arr), X)
            return arr
        end
        function $_SET!(arr1::$SETARR{N}, arr2::$SETARR{N}) where {N}
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

Create two methods to avoid method ambiguities for a lazy set operation with
respect to neutral-element and absorbing-element set types.

### Input

- `SET`  -- set type of lazy operation
- `NEUT` -- set type of neutral element
- `ABS`  -- set type of absorbing element

### Output

A quoted expression containing the function definitions.

### Notes

This macro is used internally in other macros.

### Examples

`@neutral_absorbing(MinkowskiSum, N, A)` creates the following methods as
quoted expressions:
* `MinkowskiSum(N, A) = A`
* `MinkowskiSum(A, N) = A`
"""
macro neutral_absorbing(SET, NEUT, ABS)
    return quote
        function $SET(::$NEUT{N}, Y::$ABS{N}) where {N}
            return Y
        end
        function $SET(Y::$ABS{N}, ::$NEUT{N}) where {N}
            return Y
        end
    end
end

"""
    @array_neutral(FUN, NEUT, SETARR)

Create two methods to avoid method ambiguities for a lazy set operation with
respect to the neutral-element set type and the array set type.

### Input

- `FUN`     -- function name
- `NEUT`    -- set type of neutral element
- `SETARR`  -- set type of array version

### Output

A quoted expression containing the function definitions.

### Examples

`@array_neutral(MinkowskiSum, N, ARR)` creates the following methods as
quoted expressions:
* `MinkowskiSum(N, ARR) = ARR`
* `MinkowskiSum(ARR, N) = ARR`
"""
macro array_neutral(FUN, NEUT, SETARR)
    return quote
        function $FUN(::$NEUT{N}, X::$SETARR{N}) where {N}
            return X
        end
        function $FUN(X::$SETARR{N}, ::$NEUT{N}) where {N}
            return X
        end
    end
end

"""
    @array_absorbing(FUN, ABS, SETARR)

Create two methods to avoid method ambiguities for a lazy set operation with
respect to the absorbing-element set type and the array set type.

### Input

- `FUN`     -- function name
- `ABS`     -- set type of absorbing element
- `SETARR`  -- set type of array version

### Output

A quoted expression containing the function definitions.

### Examples

`@array_absorbing(MinkowskiSum, ABS, ARR)` creates the following methods as
quoted expressions:
* `MinkowskiSum(ABS, ARR) = ABS`
* `MinkowskiSum(ARR, ABS) = ABS`
"""
macro array_absorbing(FUN, ABS, SETARR)
    return quote
        function $FUN(Y::$ABS{N}, ::$SETARR{N}) where {N}
            return Y
        end
        function $FUN(::$SETARR{N}, Y::$ABS{N}) where {N}
            return Y
        end
    end
end
