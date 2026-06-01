export ExactSum, ⊞

"""
    ExactSum{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}

Type that represents the exact sum of two sets [KochdumperA21; Proposition 10](@citet)

### Fields

- `X` -- set
- `Y` -- set

### Notes

The convenience aliases `⊞` is also available. `⊞` can be typed by `\\boxplus<tab>`.
"""
struct ExactSum{N,S1,S2} <: LazySet{N}
    X::S1
    Y::S2

    function ExactSum(X::S1, Y::S2) where {N,S1<:LazySet{N},S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "The sets must have the same ambient dimension."
        return new{N,S1,S2}(X, Y)
    end
end

function dim(ES::ExactSum)
    return dim(ES.X)
end

⊞(X::LazySet, Y::LazySet) = ExactSum(X, Y)

isoperationtype(::Type{<:ExactSum}) = true
concrete_function(::Type{<:ExactSum}) = exact_sum

# interface for binary set operations
first(ES::ExactSum) = ES.X
second(ES::ExactSum) = ES.Y
@declare_binary_operation(ExactSum)
