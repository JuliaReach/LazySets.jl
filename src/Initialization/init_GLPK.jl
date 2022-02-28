# define constants whose name changed with different versions
mf = Pkg.Operations.Context().env.manifest
ver = mf[findfirst(v -> v.name == "GLPK", mf)].version
if ver >= v"1"
    const GLPK_ON = GLPK.GLP_ON
else
    const GLPK_ON = GLPK.ON
end
