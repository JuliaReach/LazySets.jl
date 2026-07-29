module SetProgExt

using JuMP: Model, optimize!, @constraint, @objective, @variable
using LazySets: AbstractPolyhedron, chebyshev_center_radius, default_sdp_solver,
                polyhedron
using LazySets.EllipsoidModule: Ellipsoid
import SetProg
using SetProg: InteriorPoint, nth_root, value, volume
using SetProg.Sets: Translation, ellipsoid
import SetProg.Sets as Sets
import Base: convert
import LazySets.Approximations: _underapproximate_ellipsoid

function _underapproximate_ellipsoid(P::AbstractPolyhedron{N}, ::Type{Ellipsoid};
                                     backend=default_sdp_solver(),
                                     interior_point::AbstractVector{N}=chebyshev_center_radius(P)[1]) where {N}
    # create SDP model
    model = Model(backend)

    # create symbolic ellipsoid S
    point = InteriorPoint(interior_point)
    @variable(model, S, SetProg.Ellipsoid(point=point))

    # add subset constraint
    Q = polyhedron(P)  # convert P to a Polyhedra.polyhedron
    @constraint(model, S ⊆ Q)

    # add maximum-volume objective
    @objective(model, Max, nth_root(volume(S)))

    # solve SDP
    optimize!(model)

    # obtain solution
    set = value(S)

    # convert to SetProg representation
    E = ellipsoid(set)

    # convert result to a LazySets.Ellipsoid
    return convert(Ellipsoid, E)
end

function convert(::Type{<:Ellipsoid}, TE::Translation{<:Sets.Ellipsoid})
    return Ellipsoid(TE.c, inv(TE.set.Q))
end

function convert(::Type{<:Ellipsoid}, E::Sets.Ellipsoid)
    return Ellipsoid(inv(E.Q))
end

end  # module
