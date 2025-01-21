# https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane#Closest_point_and_distance_for_a_hyperplane_and_arbitrary_point
@commutative function distance(x::AbstractVector, H::Hyperplane; p::Real=2)
    @assert length(x) == dim(H) "incompatible dimensions $(length(x)) and $(dim(H))"

    # wrong!
    return norm(abs(dot(x, H.a) - H.b) / dot(H.a, H.a) * H.a, p)
end
