export sign_cadlag

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
    @neutral(SET, NEUT)

Creates functions to make a set type behave commutative with a given neutral
element set type.

### Input

- `SET` -- set type

### Output

Nothing.

### Notes

This macro generates three functions (and possibly two more if `@absorbing` has
been used in advance).

### Examples

`@neutral(MinkowskiSum, N)` creates the following functions:
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
        if isdefined(:absorbing) && method_exists(absorbing, (Type{$SET},))
            ABS = absorbing($SET)
            function $SET(::$NEUT{N}, Y::ABS{N}) where {N<:Real}
                return Y
            end
            function $SET(Y::ABS{N}, ::$NEUT{N}) where {N<:Real}
                return Y
            end
        end
    end
    return nothing
end

"""
    @absorbing(SET, ABS)

Creates functions to make a set type behave commutative with a given absorbing
element set type.

### Input

- `SET` -- set type

### Output

Nothing.

### Notes

This macro generates three functions (and possibly two more if `@absorbing` has
been used in advance).

### Examples

`@absorbing(MinkowskiSum, A)` creates the following functions:
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
        if isdefined(:neutral) && method_exists(neutral, (Type{$SET},))
            NEUT = neutral($SET)
            function $SET(::NEUT{N}, Y::$ABS{N}) where {N<:Real}
                return Y
            end
            function $SET(Y::$ABS{N}, ::NEUT{N}) where {N<:Real}
                return Y
            end
        end
    end
    return nothing
end

# TODO document
#
# TODO add ambiguity functions also to the other macros such that the order of
# calling the macros does not matter
#
# TODO remove now redundant function from the three set files
macro declare_array_version(SET, SETARR)
    @eval begin
        # create function to obtain the array version
        function array_constructor(::Type{$SET})
            return $SETARR
        end

        # create functions to use the array version functions
        function $SET(X::LazySet{N}, arr::$SETARR{N}) where {N<:Real}
            return $SETARR(arr, X)
        end
        function $SET(arr::$SETARR{N}, X::LazySet{N}) where {N<:Real}
            return $SETARR(arr, X)
        end
        function $SET(arr1::$SETARR{N}, arr2::$SETARR{N}) where {N<:Real}
            return $SETARR(arr1, arr2)
        end

        # create functions for array version
        function $SETARR(X::LazySet{N}, arr::$SETARR{N}) where {N<:Real}
            push!(array(arr), X)
            return arr
        end
        function $SETARR(arr::$SETARR{N}, X::LazySet{N}) where {N<:Real}
            push!(array(arr), X)
            return arr
        end
        function $SETARR(arr1::$SETARR{N}, arr2::$SETARR{N}) where {N<:Real}
            append!(array(arr1), array(arr2))
            return arr1
        end

        # handle method ambiguities with neutral elements
        if isdefined(:neutral) && method_exists(neutral, (Type{$SET},))
            NEUT = neutral($SET)
            function $SET(::NEUT{N}, X::$SETARR{N}) where {N<:Real}
                return X
            end
            function $SET(X::$SETARR{N}, ::NEUT{N}) where {N<:Real}
                return X
            end
        end
        if isdefined(:neutral) && method_exists(neutral, (Type{$SETARR},))
            NEUT = neutral($SETARR)
            function $SETARR(::NEUT{N}, X::$SETARR{N}) where {N<:Real}
                return X
            end
            function $SETARR(X::$SETARR{N}, ::NEUT{N}) where {N<:Real}
                return X
            end
        end
        # handle method ambiguities with absorbing elements
        if isdefined(:absorbing) && method_exists(absorbing, (Type{$SET},))
            ABS = absorbing($SET)
            function $SET(Y::ABS{N}, ::$SETARR{N}) where {N<:Real}
                return Y
            end
            function $SET(::$SETARR{N}, Y::ABS{N}) where {N<:Real}
                return Y
            end
        end
        if isdefined(:absorbing) && method_exists(absorbing, (Type{$SETARR},))
            ABS = absorbing($SETARR)
            function $SETARR(Y::ABS{N}, ::$SETARR{N}) where {N<:Real}
                return Y
            end
            function $SETARR(::$SETARR{N}, Y::ABS{N}) where {N<:Real}
                return Y
            end
        end
    end
end
