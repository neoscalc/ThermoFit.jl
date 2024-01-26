using ThermoFit

CST = global_constants()

 # Test bingo (MATLAB result = 69.277480978266283)
obs_comp = [1,1.2,3.1]
obs_unc = [0.1,0.023,0.3]
mod_comp = [0.9,1.3,3.4]

bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)

# Test generate uncertainties from Bingo
obs_comp = [1,1.2,3.1]
obs_unc = bingo_generate_fake_uncertainties(obs_comp)
mod_comp = [0.9,1.3,3.4]

bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)

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



# test MAGEMin and Bingo with a manually defined biotite composition -----------------------------------------------------------
using ThermoFit
CST = global_constants()
PARAMS = global_parameters()

out = callMAGEMin()

# Stop if biotite is not in the list of phases
if !("bi" in out.ph)
    println("biotite not in the list of phases")
    return
end

comp_structural_formula_clean, oxides = calc_structural_formula_element_from_output(out,"bi",12)

constraint_composition = [2.60, 1.60, 1.2, 1.6, 0.98, 0.12,0.01]
constraint_element = ["Si","Al","Mg", "Fe", "K", "Ti", "Mn"]

comp_structural_formula_clean_ordered = fix_order_structural_formula(comp_structural_formula_clean, oxides, constraint_element)

# Generate fake uncertainties
constraint_uncertainties = bingo_generate_fake_uncertainties(constraint_composition)

# Call bingo_calculate_qcmp_phase
qcmp_phase = bingo_calculate_qcmp_phase(comp_structural_formula_clean_ordered,constraint_composition,constraint_uncertainties)




# test objective_function ---------------------------------------------------------------------------------------------------------
using ThermoFit
CST = global_constants()
PARAMS = global_parameters()

database = "mp"

constraint_A = Constraint(8.0, 650, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])
constraint_B = Constraint(8.0, 600, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])

constraints = [constraint_A, constraint_B]

phase = "bi"
margules = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

@btime objective_function(margules, constraints, database, phase, "wt")
# 213.704 ms (896 allocations: 58.55 KiB)








# test inversion_run---------------------------------------------------------------------------------------------------------

using ThermoFit

w_names =  ["W(phl,annm)",
            "W(phl,obi)",
            "W(phl,east)",
            "W(phl,tbi)",
            "W(phl,fbi)",
            "W(phl,mmbi)",
            "W(annm,obi)",
            "W(annm,east)",
            "W(annm,tbi)",
            "W(annm,fbi)",
            "W(annm,mmbi)",
            "W(obi,east)",
            "W(obi,tbi)",
            "W(obi,fbi)",
            "W(obi,mmbi)",
            "W(east,tbi)",
            "W(east,fbi)",
            "W(east,mmbi)",
            "W(tbi,fbi)",
            "W(tbi,mmbi)",
            "W(fbi,mmbi)"]

w_initial_values = [12  0  0 ; 
                    4  0  0 ; 
                    10  0  0 ; 
                    30  0  0 ; 
                    8  0  0 ; 
                    9  0  0 ; 
                    8  0  0 ; 
                    15  0  0 ; 
                    32  0  0 ; 
                    13.6  0  0 ; 
                    6.3  0  0 ; 
                    7  0  0 ; 
                    24  0  0 ; 
                    5.6  0  0 ; 
                    8.1  0  0 ; 
                    40  0  0 ; 
                    1  0  0 ; 
                    13  0  0 ; 
                    40  0  0 ; 
                    30  0  0 ; 
                    11.6  0  0]

w_lower_bounds =   [0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0]

w_upper_bounds =   [0 0 0;
                    0 0 0;
                    60 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0]

thermodynamic_database = "mp"
solid_solution = "bi"

# Create the JOB structure
JOB = job(thermodynamic_database, solid_solution, w_names, w_initial_values, w_lower_bounds, w_upper_bounds)

# Check JOB
job_check_consistency(JOB)

# Import the constraints
CST = global_constants()
PARAMS = global_parameters()

path_mineral = "test/data/biotite.csv"
path_bulk = "test/data/bulk.csv"
path_pt = "test/data/pt.csv"

constraints = load_constraints(path_bulk, path_mineral, path_pt,["Si","Al","Mg", "Fe", "K", "Ti", "Mn"])

# Call the inversion subroutine 
res, norm = inversion_run(JOB, constraints)

# println(JOB.w_initial_values) # Note that we update JOB.w_initial_values in this version
println(res.minimizer.*norm)

# 1

