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
