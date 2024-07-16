function ρ(d::AbstractVector, P::Polygon)
    require(@__MODULE__, :LazySets; fun_name="ρ")

    return _ρ_vertices(d, P.vertices)
end
