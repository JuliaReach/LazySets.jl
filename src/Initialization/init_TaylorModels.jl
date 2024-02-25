using .TaylorModels: domain
using .TaylorModels.TaylorSeries: numtype

_eltype_TM(TMi) = numtype(polynomial(TMi))

# check that a vector of Taylor models has the [-1, 1] domain
function _has_normalized_domain(vTM::Vector)
    return all(TMi -> all(==(sym_itv(_eltype_TM(TMi))), domain(TMi)), vTM)
end

eval(load_taylormodels_convert_polynomial_zonotope())
