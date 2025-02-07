module ThermoFit
using Base.Threads
using Reexport: @reexport
using Parameters
using Optim
using ProgressBars
using OrderedCollections
using YAML

@reexport using MAGEMin_C

include("constraint.jl")
include("params.jl")
include("bingo.jl")
include("forward.jl")
include("inversion.jl")
include("resplot.jl")

export Constraint, Constraint_legacy, read_constraints_from_yml

export global_constants, global_parameters, JOB
export utils, load_bulk, calc_structural_formula_element_from_output, fix_order_structural_formula, load_mineral_composition, load_bulk_composition, load_constraints, load_pt, job_check_consistency, get_variables_optim, calculate_w_g
export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export forward, callMAGEMin
export JOB, variable_optimised, inversion_run, objective_function
export resplot, pixelmap_calc

end