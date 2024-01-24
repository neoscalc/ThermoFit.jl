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
using ThermoFit

# CST = global_constants()

database = "mp";                    # select database here, ig, igd, alk, mp, mb, um

global gv, z_b, DB, splx_data   = init_MAGEMin(database);

gv = use_predefined_bulk_rock(gv, 0, database);

# Print information on the selected database
gv, z_b, DB, splx_data = pwm_init(5, 650, gv, z_b, DB, splx_data);
ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

print("\n---------------------------------\n");
print("    Database information (",database,")\n");
print("----------------------------------\n");
for i=1:gv.len_ss
    print("   ",ss_names[i],": ",i,"; n_W's, ",ss_struct[i].n_w,"; n_em's, ",ss_struct[i].n_em,"\n")
end
print("----------------------------------\n\n");

# Test 1: run MAGEMin with default parameters
gv, z_b, DB, splx_data = pwm_init(5, 650, gv, z_b, DB, splx_data);
out       = pwm_run(gv, z_b, DB, splx_data);


# Test 2: run MAGEMin and display the Margules for biotite
gv, z_b, DB, splx_data = pwm_init(5, 650, gv, z_b, DB, splx_data);

# get the solution phase structure (size gv.len_ss)
ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
ss = 3
W  = unsafe_wrap(Vector{Cdouble},ss_struct[ss].W, ss_struct[ss].n_w)
#out = pwm_run(gv, z_b, DB, splx_data);

# Test 3: run MAGEMin and modify the Margules for biotite
gv, z_b, DB, splx_data = pwm_init(5, 650, gv, z_b, DB, splx_data);
ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
ss = 3
W  = unsafe_wrap(Vector{Cdouble},ss_struct[ss].W, ss_struct[ss].n_w);
W[2]  = 30
out = pwm_run(gv, z_b, DB, splx_data);

# out.SS_vec[findfirst(x->x=="bi", out.ph)].Comp ./sum(out.SS_vec[findfirst(x->x=="bi", out.ph)].Comp)  .* 22
finalize_MAGEMin(gv,DB, z_b)


using ThermoFit
CST = global_constants()
PARAMS = global_params()

out = callMAGEMin()

comp_structural_formula_clean = calc_structural_formula_element_from_output(out,"bi",12)

constraint_composition = [2.73, 1.27, 0, 1.2, 1.6, 0.98, 0, 0.12, 0.01]

# Generate fake uncertainties
constraint_uncertainties = bingo_generate_fake_uncertainties(constraint_composition)

# Call bingo_calculate_qcmp_phase
qcmp_phase = bingo_calculate_qcmp_phase(comp_structural_formula_clean,constraint_composition,constraint_uncertainties)




# Create an empty vector of Int64 with length n_em
# 

















# function calc_structural_formula_element_from_output(out,ss_name,oxygen_norm)
    
#     # find the indices of the oxides in CTS.oxides_definition (Philip)
#     oxide_idx = zeros(length(out.oxides))
#     for i = 1:length(out.oxides)
#         oxide_idx[i] = CST.oxide_index[out.oxides[i]]
#     end

#     println(oxide_idx)


#     ss_idx = findfirst(x->x=="bi", out.ph)
#     oxide_mol_comp = out.SS_vec[ss_idx].Comp

#     # println(oxide_mol_comp)


# end



# function meth1(out,ss_name,oxygen_norm)
#     # find the indices of the oxides in CTS.oxides_definition (Pierre)
#     oxide_idx = [findfirst(x->x==oxide, CST.oxides_definition) for oxide in out.oxides]
#     println(oxide_idx)
# end
# function meth2(out,ss_name,oxygen_norm)
#     # find the indices of the oxides in CTS.oxides_definition (Philip)
#     oxide_idx_2 = zeros(length(out.oxides))
#     for i = 1:length(out.oxides)
#         oxide_idx_2[i] = CST.oxide_index[out.oxides[i]]
#     end
#     println(oxide_idx_2)
# end

# @time meth1(out,"bi",12)
# @time meth2(out,"bi",12)