using .SparsePolynomialZonotopeModule: SparsePolynomialZonotope, exact_sum

export ExactSum, ⊞

"""
    ExactSum{N, S1<:SparsePolynomialZonotope{N}, S2<:SparsePolynomialZonotope{N}} <: LazySet{N}

Type that represents the exact sum of two sparse polynomial zonotopes [KochdumperA21; Proposition 10](@citet)

### Fields

- `X` -- set
- `Y` -- set

### Notes

The convenience aliases `⊞` is also available. `⊞` can be typed by `\\boxplus<tab>`.
"""
struct ExactSum{N,S1,S2} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function ExactSum(X::S1,
                      Y::S2) where {N,
                                    S1<:SparsePolynomialZonotope{N},
                                    S2<:SparsePolynomialZonotope{N}}
        @assert dim(X) == dim(Y)
        return ExactSum{N,S1,S2}(X, Y)
    end
end

⊞(X::SparsePolynomialZonotope, Y::SparsePolynomialZonotope) = ExactSum(X, Y)

isoperationtype(::Type{<:ExactSum}) = true
concrete_function(::Type{<:ExactSum}) = exact_sum

function concretize(ES::ExactSum)
    return exact_sum(ES.X, ES.Y)
end
