using ThermoFit

CST = global_constants()

 # Test bingo (MATLAB result = 69.277480978266283)
obs_comp = [1,1.2,3.1]
obs_unc = [0.1,0.023,0.3]
mod_comp = [0.9,1.3,3.4]

bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc,CST)

# Test generate uncertainties from Bingo
obs_comp = [1,1.2,3.1]
obs_unc = bingo_generate_fake_uncertainties(obs_comp)
mod_comp = [0.9,1.3,3.4]

bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc,CST)

# Test load_bulk
cd(@__DIR__)
bulk = load_bulk("data/bulk.csv")
print(bulk[1, :])


# Test call MAGEMin
using MAGEMin_C

database = "mp";                    # select database here, ig, igd, alk, mp, mb, um

global gv, z_b, DB, splx_data   = init_MAGEMin(database)

# get the solution phase structure (size gv.len_ss)
# ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

# ss = 6

# W  = unsafe_wrap(Vector{Cdouble},ss_struct[ss].W, ss_struct[ss].n_w)

out       = pwm_run(gv, z_b, DB, splx_data);
