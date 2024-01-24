module ThermoFit

using Reexport
using Parameters
using CSV
using DataFrames
using MAGEMin_jll
@reexport using MAGEMin_C

include("bingo.jl")
include("forward.jl")
include("inversion.jl")
include("params.jl")
include("utils.jl")





export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export forward, callMAGEMin
export Constraint, objective_function
export global_constants, global_parameters
export load_bulk
export utils




end