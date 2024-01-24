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
using MAGEMin_C

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

out.SS_vec[findfirst(x->x=="bi", out.ph)].Comp ./sum(out.SS_vec[findfirst(x->x=="bi", out.ph)].Comp)  .* 22

calc_structural_formula_element_from_output(out,"bi",12)

finalize_MAGEMin(gv,DB, z_b)



function calc_structural_formula_element_from_output(out,ss_name,oxygen_norm)

    oxides = out.oxides
    ss_idx = findfirst(x->x=="bi", out.ph)

    oxide_mol_comp = out.SS_vec[ss_idx].Comp

    println(oxide_mol_comp)


end





test_tuple::Tuple{Vector{String},Vector{Int64},Vector{Int64}}= ("test", 1, 1)