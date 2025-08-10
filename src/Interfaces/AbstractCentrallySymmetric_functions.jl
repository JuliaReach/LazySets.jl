# common functions of AbstractCentrallySymmetric, which are also shared with
# AbstractCentrallySymmetricPolytope because Julia does not have multiple inheritance

const ACS = Union{AbstractCentrallySymmetric,AbstractCentrallySymmetricPolytope}

@inline function dim(S::ACS)
    return length(center(S))
end

for T in Base.uniontypes(ACS)
    @eval begin
        """
        # Extended help

            an_element(S::$($T))

        ### Output

        The center of the centrally symmetric set.
        """
        function an_element(S::$T)
            return center(S)
        end
    end
end

function isempty(X::ACS, witness::Bool=false)
    return witness ? (false, an_element(X)) : false
end

for T in Base.uniontypes(ACS)
    @eval begin
        """
        # Extended help

            isuniversal(S::$($T), [witness]::Bool=false)

        ### Algorithm

        A witness is obtained by computing the support vector in direction
        `d = [1, 0, …, 0]` and adding `d` on top.
        """
        function isuniversal(S::$T, witness::Bool=false)
            if witness
                N = eltype(S)
                d = SingleEntryVector{N}(1, dim(S))
                w = σ(d, S) + d
                return (false, w)
            else
                return false
            end
        end
    end
end

@inline function center(S::ACS, i::Int)
    return center(S)[i]
end

for T in Base.uniontypes(ACS)
    @eval begin
        """
        # Extended help

            extrema(S::$($T))

        ### Notes

        The result is equivalent to `(low(S), high(S))`.

        ### Algorithm

        We compute `high(S)` and then compute the lowest coordinates with the help of
        `center(S)` (which is assumed to be cheaper to obtain).
        """
        function extrema(S::$T)
            # h = c + r
            h = high(S)
            # l = c - r = -c - r + 2 * c = 2 * c - h
            l = 2 .* center(S) .- h
            return (l, h)
        end
    end
end

for T in Base.uniontypes(ACS)
    @eval begin
        """
        # Extended help

            extrema(S::$($T), i::Int)

        ### Notes

        The result is equivalent to `(low(S, i), high(S, i))`.

        ### Algorithm

        We compute `high(S, i)` and then compute the lowest coordinates with the help of
        `center(S, i)` (which is assumed to be cheaper to obtain).
        """
        function extrema(S::$T, i::Int)
            # h = c + r
            h = high(S, i)
            # l = c - r = -c - r + 2 * c = 2 * c - h
            l = 2 * center(S, i) - h
            return (l, h)
        end
    end
end
