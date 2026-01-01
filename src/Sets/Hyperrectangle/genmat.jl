function genmat(H::Hyperrectangle{N,<:AbstractVector,<:SparseVector{N}}) where {N}
    n = dim(H)
    nze, nzv = findnz(H.radius)
    return sparse(nze, 1:length(nze), nzv, n, length(nze))
end

function load_StaticArraysCore_genmat_Hyperrectangle()
    return quote
        using .StaticArraysCore: SVector

        function genmat(H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {L,N}
            gens = zeros(StaticArraysCore.MMatrix{L,L,N})
            j = 0
            @inbounds for i in 1:L
                ri = radius_hyperrectangle(H, i)
                if !isapproxzero(ri)
                    j += 1
                    gens[i, j] = ri
                end
            end
            return StaticArraysCore.SMatrix{L,j}(view(gens, :, 1:j))
        end

        # this function is type stable, but it does not prune the generators
        # according to flat dimensions of H
        function _genmat_static(H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {L,N}
            gens = zeros(StaticArraysCore.MMatrix{L,L,N})
            @inbounds for i in 1:L
                ri = radius_hyperrectangle(H, i)
                gens[i, i] = ri
            end
            return StaticArraysCore.SMatrix{L,L}(gens)
        end
    end
end  # load_StaticArraysCore_genmat_Hyperrectangle
