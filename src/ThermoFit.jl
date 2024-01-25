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


export global_constants, global_parameters, job
export utils, load_bulk, calc_structural_formula_element_from_output, fix_order_structural_formula, load_mineral_composition, load_bulk_composition, load_constraints, load_pt, job_check_consistency, get_variables_optim, calculate_w_g
export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export forward, callMAGEMin
export inversion, constraint, objective_function

end