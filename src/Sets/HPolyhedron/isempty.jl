# this method is required mainly for HPolytope (because the fallback for
# AbstractPolytope is incorrect with no constraints)
#
# the method also treats a corner case for problems with Rationals in LP solver
function isempty(P::HPoly,
                 witness::Bool=false;
                 use_polyhedra_interface::Bool=false,
                 solver=nothing,
                 backend=nothing)
    if length(constraints_list(P)) < 2
        return witness ? (false, an_element(P)) : false
    end
    return _isempty_polyhedron(P, witness;
                               use_polyhedra_interface=use_polyhedra_interface,
                               solver=solver, backend=backend)
end
