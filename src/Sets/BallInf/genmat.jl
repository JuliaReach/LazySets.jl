function load_genmat_ballinf_static()
    return quote
        function genmat(B::BallInf{N,SVector{L,N}}) where {L,N}
            if isflat(B)
                return SMatrix{L,0,N,0}()
            else
                gens = zeros(MMatrix{L,L})
                @inbounds for i in 1:L
                    gens[i, i] = B.radius
                end
                return SMatrix(gens)
            end
        end
    end
end  # quote / load_genmat_ballinf_static()
