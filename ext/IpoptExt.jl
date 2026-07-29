module IpoptExt

# TODO undo when dropping that version (works in v"1.12" at least; not tested in v"1.11")
# using Ipopt: Optimizer  # inlined below because of a loading issue in v"1.10"
import Ipopt
import LazySets: _default_nln_solver

function _default_nln_solver(N::Type{<:Real})
    return Ipopt.Optimizer
end

end  # module
