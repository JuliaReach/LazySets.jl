# default semidefinite-programming solver
function default_sdp_solver(N::Type{<:Number})
    error("no default SDP solver for numeric type $N")
end
