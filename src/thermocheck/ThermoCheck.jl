module ThermoCheck

using ..ThermoFit
using CairoMakie
using ProgressBars

include("pixelmap.jl")
include("qfactor.jl")

export pixelmap, pixmap_misfit_Ti_in_Bt_temperatures, pixmap_misfit_Ti_in_Bt_temperatures, pixmap_phase_comp_og_v_mod

end