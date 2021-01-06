export cartesian_product

"""
    cartesian_product(P1::VPolytope, P2::VPolytope;
                      [backend]=default_polyhedra_backend(P1))

Compute the Cartesian product of two polytopes in V-representation.

### Input

- `P1`      -- polytope
- `P2`      -- another polytope
- `backend` -- (optional, default: `default_polyhedra_backend(P1)`) the
               backend for polyhedral computations; see [Polyhedra's
               documentation](https://juliapolyhedra.github.io/) for further
               information

### Output

The `VPolytope` obtained by the concrete Cartesian product of `P1` and `P2`.
"""
function cartesian_product(P1::VPolytope, P2::VPolytope;
                           backend=default_polyhedra_backend(P1))
    require(:Polyhedra; fun_name="cartesian_product")
    Pcp = Polyhedra.vcartesianproduct(polyhedron(P1; backend=backend),
                                      polyhedron(P2; backend=backend))
    return VPolytope(Pcp)
end

"""
    cartesian_product(P1::HPoly, P2::HPoly;
                      [backend]=default_polyhedra_backend(P1)
                     )

Compute the Cartesian product of two polyhedra in H-representaion.

### Input

- `P1`         -- polyhedron
- `P2`         -- another polyhedron
- `backend`    -- (optional, default: `default_polyhedra_backend(P1)`)
                  the polyhedral computations backend

### Output

The polyhedron obtained by the concrete cartesian product of `P1` and `P2`.

### Notes

For further information on the supported backends see
[Polyhedra's documentation](https://juliapolyhedra.github.io/).
"""
function cartesian_product(P1::HPoly, P2::HPoly;
                           backend=default_polyhedra_backend(P1)
                          )
    require(:Polyhedra; fun_name="`cartesian_product")
    Pcp = Polyhedra.hcartesianproduct(polyhedron(P1; backend=backend),
                                      polyhedron(P2; backend=backend))
    return convert(basetype(P1), Pcp)
end
