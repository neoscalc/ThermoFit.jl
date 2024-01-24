module ThermoFit

using Reexport
using Parameters
using CSV
using DataFrames
using MAGEMin_jll
@reexport using MAGEMin_C

include("params.jl")
include("bingo.jl")
include("forward.jl")
include("inversion.jl")
include("utils.jl")





export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export forward, callMAGEMin
export Constraint, objective_function
export global_constants, global_params
export load_bulk
export utils, calc_structural_formula_element_from_output




end