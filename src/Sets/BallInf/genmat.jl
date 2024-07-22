function load_StaticArraysCore_genmat()
    return quote
        using StaticArraysCore: SVector

        function genmat(B::BallInf{N,SVector{L,N}}) where {L,N}
            if isflat(B)
                return StaticArraysCore.SMatrix{L,0,N,0}()
            else
                gens = zeros(StaticArraysCore.MMatrix{L,L})
                @inbounds for i in 1:L
                    gens[i, i] = B.radius
                end
                return StaticArraysCore.SMatrix(gens)
            end
        end
    end
end  # load_StaticArraysCore_genmat
