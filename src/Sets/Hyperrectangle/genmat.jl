function genmat(H::Hyperrectangle{N,<:AbstractVector,<:SparseVector{N}}) where {N}
    n = dim(H)
    nze, nzv = findnz(H.radius)
    return sparse(nze, 1:length(nze), nzv, n, length(nze))
end

function load_StaticArraysCore_genmat()
    return quote
        using .StaticArraysCore: SVector

        function genmat(H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {L,N}
            gens = zeros(StaticArraysCore.MMatrix{L,L,N})
            nzcol = Vector{Int}()
            @inbounds for i in 1:L
                r = H.radius[i]
                if !isapproxzero(r)
                    gens[i, i] = r
                    push!(nzcol, i)
                end
            end
            m = length(nzcol)
            return StaticArraysCore.SMatrix{L,m}(view(gens, :, nzcol))
        end

        # this function is type stable, but it does not prune the generators
        # according to flat dimensions of H
        function _genmat_static(H::Hyperrectangle{N,SVector{L,N},SVector{L,N}}) where {L,N}
            gens = zeros(StaticArraysCore.MMatrix{L,L,N})
            @inbounds for i in 1:L
                r = H.radius[i]
                gens[i, i] = r
            end
            return StaticArraysCore.SMatrix{L,L}(gens)
        end
    end
end  # load_StaticArraysCore_genmat
