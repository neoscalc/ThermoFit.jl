module ThermoFit

using Reexport: @reexport
using Random
using Dates
using Statistics
using Optim
using ProgressBars
using OrderedCollections
using CSV
using DataFrames
using YAML
using CairoMakie

@reexport using Base.Threads
@reexport using MAGEMin_C

include("constraint.jl")
include("forward.jl")
include("inversion.jl")
include("logging.jl")
include("loss.jl")
include("objective_f.jl")
include("plot.jl")

include("thermocheck/ThermoCheck.jl")
using .ThermoCheck
export ThermoCheck

export Constraint, read_constraints_from_yml, gen_constraints_for_functional_inv, read_FPWMP_bulks
export forward_call, calculate_w_g
export JOB, print_job, variable_optimised, inversion
export print_job, print_constraints, print_results, find_optimal_solution
export objective_function, objective_function_func_relation
export loss_Qfactor, quality_factor, bingo_generate_fake_uncertainties, chi_squared, Ti_in_Bt_misfit, Ti_sat_misfit
export plot_convergence, x_optim_v_resiudal
export pixelmap, pixmap_misfit_Ti_in_Bt_temperatures, pixmap_misfit_Ti_in_Bt_temperatures, pixmap_phase_comp_og_v_mod

end