@validate function project(L::Line{N}, block::AbstractVector{Int}; kwargs...) where {N}
    d = L.d[block]
    if iszero(d)
        require(@__MODULE__, :LazySets; fun_name="project")

        return Singleton(L.p[block])  # projected out all nontrivial dimensions
    elseif length(d) == 1
        require(@__MODULE__, :LazySets; fun_name="project")

        return Universe{N}(1)  # special case: 1D line is a universe
    else
        return Line(L.p[block], d)
    end
end
