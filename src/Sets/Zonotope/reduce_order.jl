function load_StaticArraysCore_reduce_order()
    return quote
        using .StaticArraysCore: MMatrix
        using ..LazySets: AbstractReductionMethod

        # conversion for static matrix
        function reduce_order(Z::Zonotope{N,<:AbstractVector,<:MMatrix}, r::Real,
                              method::AbstractReductionMethod=GIR05()) where {N}
            return reduce_order(Zonotope(center(Z), StaticArraysCore.SMatrix(genmat(Z))), r, method)
        end
    end
end # load_StaticArraysCore_reduce_order
