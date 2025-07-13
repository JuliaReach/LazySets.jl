@validate function σ(d::AbstractVector, P::Polygon)
    require(@__MODULE__, :LazySets; fun_name="σ")

    return _σ_vertices(d, P.vertices)
end
