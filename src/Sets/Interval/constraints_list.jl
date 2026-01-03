function constraints_list(X::Interval)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    N = eltype(X)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2)
    e₁ = SingleEntryVector(1, 1, one(N))
    @inbounds constraints[1] = HalfSpace(e₁, max(X))
    @inbounds constraints[2] = HalfSpace(-e₁, -min(X))
    return constraints
end

function _constraints_list_Vector(X::Interval)
    require(@__MODULE__, :LazySets; fun_name="constraints_list")

    N = eltype(X)
    constraints = Vector{HalfSpace{N,Vector{N}}}(undef, 2)
    @inbounds constraints[1] = HalfSpace([one(N)], max(X))
    @inbounds constraints[2] = HalfSpace([-one(N)], -min(X))
    return constraints
end
