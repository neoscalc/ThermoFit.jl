module ThermoFit

using Parameters

include("bingo.jl")
include("forward.jl")

export global_constants
export bingo, bingo_calculate_qcmp_phase, bingo_generate_fake_uncertainties

end