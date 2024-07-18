function load_reduce_order_static_zonotope()
    return quote
        using ..LazySets: AbstractReductionMethod

        # conversion for static matrix
        function reduce_order(Z::Zonotope{N,<:AbstractVector,<:MMatrix}, r::Real,
                              method::AbstractReductionMethod=GIR05()) where {N}
            return reduce_order(Zonotope(center(Z), SMatrix(genmat(Z))), r, method)
        end
    end
end # quote / load_reduce_order_static_zonotope
