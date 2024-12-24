function reflect(P::VPolytope)
    return VPolytope(-P.vertices)
end
