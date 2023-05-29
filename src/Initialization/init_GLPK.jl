# define constants whose name changed with different versions
if isdefined(GLPK, :GLP_ON)
    const GLPK_ON = GLPK.GLP_ON
else
    const GLPK_ON = GLPK.ON
end
