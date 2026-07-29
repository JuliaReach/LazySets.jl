function default_nln_solver(N::Type{<:Real}=Float64)
    return _default_nln_solver(N)
end

# see ext/IpoptExt.jl
function _default_nln_solver(N)
    mod = Base.get_extension(@__MODULE__, :IpoptExt)
    require(mod, :Ipopt; fun_name="default_nln_solver")
    error()
end
