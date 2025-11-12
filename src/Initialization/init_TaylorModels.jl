using .TaylorModels: domain
using .TaylorModels.TaylorSeries: numtype  # NOTE: this is an internal function

_eltype_TM(TMi) = numtype(polynomial(TMi))

# check that a vector of Taylor models has the [-1, 1] domain
function _has_normalized_domain(vTM::Vector)
    return all(TMi -> all(IA.isequal_interval(sym_itv(_eltype_TM(TMi))), domain(TMi)), vTM)
end

eval(load_TaylorModels_convert_SparsePolynomialZonotope())
eval(load_TaylorModels_convert_TaylorModelN())
