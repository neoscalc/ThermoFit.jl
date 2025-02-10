using ThermoFit
using Test

@testset "ThermoFit" begin

@testset "full_inversion_test" begin
    cd(@__DIR__)
    # start with the Margules for biotite from White et al. (2014) but modify 
    # allow only the Margules W(phl,annm) to change within a range of -100 - 100 kJ/mol
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
                "W(fbi,mmbi)"];

    w_initial_values = [1.  0  0 ;
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
                        11.6  0  0];

    w_lower_bounds = copy(w_initial_values)
    w_lower_bounds[1,1] = -100
    w_upper_bounds = copy(w_initial_values)
    w_upper_bounds[1,1] = 100

    constraints = read_constraints_from_yml("test_data/gen_data_FPWMP_biotites.yml")

    # set-up two jobs for testing the two optimisation algorithms:
    thermodynamic_database = "mp";
    solid_solution = "bi";
    number_constraints_max = 10;
    number_iterations_max = 10;
    max_time_seconds = 120;

    # 1. Nelder-Mead (with normalization)
    algorithm = "NelderMead";
    normalization = true;

    job_NelderMead = JOB(thermodynamic_database, solid_solution,
                         w_names=w_names, w_initial_values=w_initial_values, w_lower_bounds=w_lower_bounds, w_upper_bounds=w_upper_bounds,
                         algorithm=algorithm, number_iterations_max=number_iterations_max, normalization=normalization,
                         number_constraints_max=number_constraints_max, max_time_seconds=max_time_seconds);
    print_job(job_NelderMead)

    # 2. Particle Swarm (without normalization)
    algorithm = "ParticleSwarm";
    normalization = false;

    job_ParticleSwarm = JOB(thermodynamic_database, solid_solution,
                            w_names=w_names, w_initial_values=w_initial_values, w_lower_bounds=w_lower_bounds, w_upper_bounds=w_upper_bounds,
                            algorithm=algorithm, number_iterations_max=number_iterations_max, normalization=normalization,
                            number_constraints_max=number_constraints_max, max_time_seconds=max_time_seconds);
    print_job(job_ParticleSwarm)

    # RUN test inversions
    res_NelderMead, norm_NelderMead = inversion(job_NelderMead, constraints)
    res_ParticleSwarm, norm_ParticleSwarm = inversion(job_ParticleSwarm, constraints)

    # check inverterted Margules parameters
    @test res_NelderMead.minimizer .* norm_NelderMead ≈ [12.]             atol=1
    # ParticleSwarm is stochastic, so we only check if the number of iterations is equal to the maximum number of iterations
    @test res_ParticleSwarm.iterations == 10
end

##############################################################################################################
# Test of individual modules and functions therein
##############################################################################################################

@testset "bingo.jl" begin

    obs_comp = [1,1.2,3.1]
    obs_unc = [0.1,0.023,0.3]
    mod_comp = [0.9,1.3,3.4]

    qcmp = bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)

    @test qcmp ≈ 69.277480978266283 atol=1e-6

end

@testset "forward.jl" begin
    phase = "bi"
    database = "mp"

    # set up the constraint
    temperature_C = 580
    pressure_GPa = 0.45
    bulk = [64.13, 0.91, 19.63, 6.85, 0.08, 2.41, 0.65, 1.38, 3.95, 40.0]
    bulk_oxides = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "H2O"]
    sys_in = "wt%"
    assemblage = ["st", "bi", "mu", "fsp", "ilm", "q", "H2O"]
    mineral_composition_apfu = Dict("bi" => [2.720008886150718, 1.559982227698564, 0.0, 0.89014494411893, 1.732872689303359, 1.0, 0.0, 0.08652624797247646, 12.0, 0.01046500475595296, 1.8269475040550471])
    mineral_elements = ["Si", "Al", "Ca", "Mg", "Fe", "K", "Na", "Ti", "O", "Mn", "H"]

    constraint = Constraint(pressure_GPa, temperature_C, bulk, bulk_oxides, sys_in, assemblage, mineral_composition_apfu, mineral_elements)

    MAGEMin_db  = Initialize_MAGEMin("mp", verbose=false)

    # identify thread and acess the MAGEMin_db of the thread
    id          = Threads.threadid()
    gv          = MAGEMin_db.gv[id]
    z_b         = MAGEMin_db.z_b[id]
    DB          = MAGEMin_db.DB[id]
    splx_data   = MAGEMin_db.splx_data[id]

    # use unaltered w_g = W_H from the mp database
    w_initial_values = [12., 4, 10, 30, 8, 9, 8, 15, 32, 13.6, 6.3, 7, 24, 5.6, 8.1, 40, 1, 13, 40, 30, 11.6]
    w_mod = [20., 4, 10, 30, 8, 9, 8, 15, 32, 13.6, 6.3, 7, 24, 5.6, 8.1, 40, 1, 13, 40, 30, 11.6]
    g0_corr_endmembers = ["phl"]
    g0_corr_initial_values = [0.0]
    g0_corr_mod = [2.0]

    # test the forward_call() function with no w_g and g0_corr
    out = forward_call(phase, database, constraint, gv, z_b, DB, splx_data; w_g = nothing, g0_corr_vec = nothing, g0_corr_em = nothing)
    @test sort(out.ph) == sort(constraint.assemblage)
    @test out.SS_vec[findfirst(out.ph .== "bi")].Comp_apfu ≈ constraint.mineral_composition_apfu["bi"] atol=1e-4

    # test the forward_call() function with unaltered w_g and g0_corr
    out = forward_call(phase, database, constraint, gv, z_b, DB, splx_data; w_g = w_initial_values, g0_corr_vec = g0_corr_initial_values, g0_corr_em = g0_corr_endmembers)
    @test sort(out.ph) == sort(constraint.assemblage)
    @test out.SS_vec[findfirst(out.ph .== "bi")].Comp_apfu ≈ constraint.mineral_composition_apfu["bi"] atol=1e-4

    # test the forward_call() function with altered w_g
    out_w_mod = forward_call(phase, database, constraint, gv, z_b, DB, splx_data; w_g = w_mod)
    @test sort(out_w_mod.ph) == sort(["mu", "bi", "st", "fsp", "g", "ilm", "q", "H2O"])
    @test out_w_mod.SS_vec[findfirst(out_w_mod.ph .== "bi")].Comp_apfu ≈ [2.703201966718807, 1.5935960665623852, 0.0, 0.8769929425165847, 1.7275929039130509, 1.0000000000000002, 0.0, 0.08873358434595208, 12.0, 0.00988253594321993, 1.822532831308096] atol=1e-4

    # test the forward_call() function with altered g0_corr
    out_g0_mod = forward_call(phase, database, constraint, gv, z_b, DB, splx_data; g0_corr_vec = g0_corr_mod, g0_corr_em = g0_corr_endmembers)
    @test sort(out_g0_mod.ph) == sort(["fsp", "bi", "mu", "st", "ilm", "g", "q", "H2O"])
    @test out_g0_mod.SS_vec[findfirst(out_g0_mod.ph .== "bi")].Comp_apfu ≈ [2.711326941147552, 1.5773461177048962, 0.0, 0.8795051695973085, 1.7337723613120772, 1.0, 0.0, 0.08775150164407976, 12.0, 0.010297908594086575, 1.8244969967118405] atol=1e-4

    #test the forward_call() function with altered w_g and g0_corr
    out_w_g0_mod = forward_call(phase, database, constraint, gv, z_b, DB, splx_data; w_g = w_mod, g0_corr_vec = g0_corr_mod, g0_corr_em = g0_corr_endmembers)
    @test sort(out_w_g0_mod.ph) == sort(["fsp", "mu", "st", "bi", "g", "ilm", "q", "H2O"])
    @test out_w_g0_mod.SS_vec[findfirst(out_w_g0_mod.ph .== "bi")].Comp_apfu ≈ [2.6981048656876174, 1.6037902686247656, 0.0, 0.8692846574817923, 1.7295717750076536, 1.0, 0.0, 0.0894953787332389, 12.0, 0.00975305446493273, 1.8210092425335223] atol=1e-4

    # test the calculate_w_g() function
    w_initial = [1 1 3;
                 1 1 4;
                 1 2 5]

    w_lower_bounds = [1 1 1;
                      1 1 1;
                      1 1 1]

    w_upper_bounds = [1 1 2;
                      1 1 2;
                      2 2 2]

    w_names = ["A", "B", "C"]

    job = JOB("mp", "bi", w_names=w_names, w_initial_values=w_initial, w_lower_bounds=w_lower_bounds, w_upper_bounds=w_upper_bounds)

    margules_optim = [10, 20, 30, 40, 50]
    P = 2.
    T = 500.

    w_matrix = copy(job.w_initial_values)   # copy needed for multi_threading, so every thread has its own copy of w_all
    w_coordinates = job.margules_optim_coord
    
    w_g = calculate_w_g(margules_optim, P, T, w_matrix, w_coordinates)

    w_mod = [1  1  30     # internally the function calculate_w_g() insert the values at the coordinates
             1  1  40
             10 20 50]

    w_g_mod = w_mod[:,1] .- T .* w_mod[:,2] .+ P .* w_mod[:,3]
    @test w_g == w_g_mod
end

@testset "inversion.jl" begin
    # test the JOB struct and its constructor function(s)
    # (1) test the JOB struct with only the required arguments
    job = JOB("mp", "bi", w_names=["W(phl,annm)"], w_initial_values=[12. 0 0], w_lower_bounds=[0. 0 0], w_upper_bounds=[15. 0 0])
    @test typeof(job) == JOB{String, Vector{String}, Matrix{Float64}, Int64, Bool, Nothing}

    # (2) test typeerror when passing wrong type to JOB constructor
    @test_throws MethodError JOB("mp", "bi", w_names=["W(phl,annm)"], w_initial_values=[12. 0 0], w_lower_bounds=[0 0 0], w_upper_bounds=[0 0 0])

    # (3) test the JOB constructor function for different cases of variables to optimise:
    #        - (I) optimise margules
    #        - (II) optimise g0 corrections
    #        - (III) optimise both margules and g0 corrections

    job_I = JOB("mp", "bi", w_names=["W(phl,annm)", "W(phl,obi)"],
                w_initial_values=[12. 0 0; 4  0  0], w_lower_bounds=[0. 0 0; 0 0 0], w_upper_bounds=[15. 0 0; 5 0 0])
    @test job_I.var_optim == [12.0, 4.0]
    @test job_I.var_optim_bounds == [0.0 15.0; 0.0 5.0]

    job_II = JOB("mp", "bi", g0_corr_endmembers = ["annm", "obi"],
                 g0_corr_initial_values = [2.0, 0.05], g0_corr_lower_bounds = [0.0, 0.0], g0_corr_upper_bounds = [5.0, 5.0])
    @test job_II.var_optim == [2.0, 0.05]
    @test job_II.var_optim_bounds == [0.0 5.0; 0.0 5.0]

    job_III = JOB("mp", "bi", w_names=["W(phl,annm)", "W(phl,obi)"],
                  w_initial_values=[12. 0 0; 4  0  0], w_lower_bounds=[0. 0 0; 0 0 0], w_upper_bounds=[15. 0 0; 5 0 0],
                  g0_corr_endmembers = ["annm", "obi"], g0_corr_initial_values = [2.0, 0.05], g0_corr_lower_bounds = [0.0, 0.0], g0_corr_upper_bounds = [5.0, 5.0])
    @test job_III.var_optim == [12.0, 4.0, 2.0, 0.05]
    @test job_III.var_optim_bounds == [0.0 15.0; 0.0 5.0; 0.0 5.0; 0.0 5.0]

    # (4) test the variable_optimised() function (used within JOB constructor)
    initial = [1 1 3;
               1 1 4;
               1 2 5]

    lower_bounds = [1 1 1;
                    1 1 1;
                    1 1 1]

    upper_bounds = [1 1 2;
                    1 1 2;
                    2 2 2]

    names = ["A", "B", "C"]

    var_opti, var_opti_bounds, var_opti_names,  var_opti_coord, n = variable_optimised(initial, lower_bounds, upper_bounds, names)

    @test var_opti == [1, 2, 3, 4, 5]

    # change the initial values at the coordinates
    initial[var_opti_coord] .= [10, 20, 30, 40, 50]

    init_mod = [1  1  30
                1  1  40
                10 20 50]

    @test initial == init_mod
end
end
