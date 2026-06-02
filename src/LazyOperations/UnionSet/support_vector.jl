"""
    σ(d::AbstractVector, cup::UnionSet; [algorithm]="support_vector")

Return a support vector of the union of two sets in a given direction.

### Input

- `d`         -- direction
- `cup`       -- union of two sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                 the support vector; if "support_vector", use the support
                 vector of each argument; if "support_function" use the support
                 function of each argument and evaluate the support vector of
                 only one of them

### Output

A support vector in the given direction.

### Algorithm

The support vector of the union of two sets ``X`` and ``Y`` can be obtained as
the vector that maximizes the support function of either ``X`` or ``Y``, i.e.,
it is sufficient to find the ``\\argmax(ρ(d, X), ρ(d, Y)])`` and evaluate its
support vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of ``X`` and ``Y`` and then compares the support function
using a dot product.

If the support function can be computed more efficiently, the alternative
implementation `algorithm="support_function"` can be used, which evaluates the
support function of each set directly and then calls only the support vector of
either ``X`` *or* ``Y``.
"""
@validate function σ(d::AbstractVector, cup::UnionSet; algorithm="support_vector")
    X, Y = cup.X, cup.Y

    if algorithm == "support_vector"
        σX, σY = σ(d, X), σ(d, Y)
        return dot(d, σX) > dot(d, σY) ? σX : σY

    elseif algorithm == "support_function"
        ρX, ρY = ρ(d, X), ρ(d, Y)
        return ρX > ρY ? σ(d, X) : σ(d, Y)

    else
        throw(ArgumentError("algorithm $algorithm unknown"))
    end
end
