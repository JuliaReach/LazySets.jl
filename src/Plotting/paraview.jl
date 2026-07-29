# see ext/WriteVTKExt.jl
function writevtk(X; file="output")
    mod = Base.get_extension(@__MODULE__, :WriteVTKExt)
    require(mod, :WriteVTK; fun_name="writevtk")
    error()
end
