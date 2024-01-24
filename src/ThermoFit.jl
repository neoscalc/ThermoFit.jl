module ThermoFit

using Reexport
using Parameters
using CSV
using DataFrames
using MAGEMin_jll
@reexport using MAGEMin_C

include("params.jl")
include("utils.jl")
include("bingo.jl")
include("forward.jl")
include("inversion.jl")


export global_constants, global_params
export utils, load_bulk, calc_structural_formula_element_from_output, fix_order_structural_formula
export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export forward, callMAGEMin
export inversion, constraint, objective_function

end