# see ext/MakieExt.jl
function plot3d(S; backend=default_polyhedra_backend(S), alpha=1.0,
                color=:blue, colormap=:viridis, colorrange=nothing,
                interpolate=false, overdraw=false, shading=true,
                transparency=true, visible=true)
    mod = Base.get_extension(@__MODULE__, :MakieExt)
    require(mod, [:Makie, :Polyhedra]; fun_name="plot3d")
    error()
end

# see ext/MakieExt.jl
function plot3d!(S; backend=default_polyhedra_backend(S), alpha=1.0,
                 color=:blue, colormap=:viridis, colorrange=nothing,
                 interpolate=false, overdraw=false, shading=true,
                 transparency=true, visible=true)
    mod = Base.get_extension(@__MODULE__, :MakieExt)
    require(mod, [:Makie, :Polyhedra]; fun_name="plot3d!")
    error()
end
