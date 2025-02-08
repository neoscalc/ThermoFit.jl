module ThermoFit
using Base.Threads
using Reexport: @reexport
using Parameters
using Optim
using ProgressBars
using OrderedCollections
using YAML

@reexport using MAGEMin_C

include("bingo.jl")
include("constraint.jl")
include("forward.jl")
include("inversion.jl")
include("resplot.jl")

export bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export Constraint, read_constraints_from_yml
export forward_call, calculate_w_g
export JOB, job_check_consistency, variable_optimised, inversion_run, objective_function
export resplot, pixelmap_calc

end