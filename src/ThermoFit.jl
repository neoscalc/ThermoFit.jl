module ThermoFit

using Parameters
using CSV
using DataFrames
# using MAGEMin_C

include("bingo.jl")
include("forward.jl")
include("params.jl")
include("utils.jl")





export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties
export global_constants, global_parameters
export load_bulk

end