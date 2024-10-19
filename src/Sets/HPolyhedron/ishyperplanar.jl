function ishyperplanar(P::HPolyhedron)
    clist = P.constraints
    m = length(clist)

    # check that the number of constraints is fine
    if m > 2
        # try to remove redundant constraints
        clist = remove_redundant_constraints(clist)
        if isempty(clist)
            # constraints are contradictory
            return false
        end
        m = length(clist)
    end
    if m != 2
        return false
    end

    # check that the two half-spaces are complementary
    return @inbounds iscomplement(clist[1], clist[2])
end
