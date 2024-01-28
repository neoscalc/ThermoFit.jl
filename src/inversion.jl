function inversion_run(JOB, constraints)

    # get the variables to variables_optim_coordinates
    variables_optim, variables_optim_bounds, variables_optim_coordinates = get_variables_optim(JOB)

    # println(variables_optim)

    # Set parameters for the inversion
    norm = variables_optim
    if JOB.normalization == true
        x0 = variables_optim ./ norm
    else
        x0 = variables_optim
    end
    nb_constraints = length(constraints)

    # Initiate MAGEMin
    MAGEMin_db  = Initialize_MAGEMin(JOB.thermodynamic_database, verbose=false);

    # # test call the objective function with the initial values (KEEP IT FOR TEST PURPOSES)
    # residual = objective_function(x0, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates)
    max_time_seconds = JOB.max_time_seconds

    # perform inversion using the Nelder-Mead algorithm
    if JOB.algorithm == "NelderMead"
        res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db), x0, NelderMead(), Optim.Options(time_limit = max_time_seconds))
    elseif JOB.algorithm == "ParticleSwarm"
        
        if JOB.normalization == true
            res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db), x0, ParticleSwarm(; lower = variables_optim_bounds[:,1] ./ norm, upper = variables_optim_bounds[:,2] ./ norm), Optim.Options(time_limit = max_time_seconds))  # n_particles = 0
        else
            println("x0 = ", x0)
            println("lower = ", [variables_optim_bounds[:,1]])
            println("upper = ", [variables_optim_bounds[:,2]])

            res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db), x0, ParticleSwarm(; lower = variables_optim_bounds[:,1], upper = variables_optim_bounds[:,2]), Optim.Options(time_limit = max_time_seconds))  # n_particles = 0
        end
    else
        res = 1e20
        println("Error: algorithm not recognised")
    end
    
    Finalize_MAGEMin(MAGEMin_db)

    return res, norm
end

"""
    objective_function(*args)

Add an objective function Doc here.
"""
function objective_function(x0, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db, sys_in = "mol")

    # define number of constraints to use for the inversion
    if isequal(JOB.number_constraints_max,0)
        nb_to_use = nb_constraints
    elseif JOB.number_constraints_max > nb_constraints
        nb_to_use = nb_constraints
    else
        nb_to_use = JOB.number_constraints_max
    end

    println("\n-> New iteration with ",nb_to_use," constraints <-")

    # Denormalise variables to optimise (Margules) for G-minimisation
    if JOB.normalization == true
        variables_optim_local = x0 .* norm
    else
        variables_optim_local = x0 
    end

    # Check if all parameters are within the bounds, if not return a very high residual:
    for i = eachindex(variables_optim_local)
        if variables_optim_local[i] < variables_optim_bounds[i,1]
            return 1e20
        elseif variables_optim_local[i] > variables_optim_bounds[i,2]
            return 1e20
        end
    end

    # Initiate vectors for residuals (to be minimised) and q_cpm (used as metric that is printed during the inversion)
    residual_i = zeros(nb_to_use)
    qcmp_all = zeros(nb_to_use)

    # # Precalculate the w_g for each constraint
    # w_g_all = zeros(nb_to_use, length(JOB.w_initial_values[:,1]))
    # for i = eachindex(nb_to_use)
    #     # Calculate w_g
    #     w_g_all[i, :] = calculate_w_g(variables_optim_local, variables_optim_coordinates, constraints[i].pressure, constraints[i].temperature, JOB)
    # end

    @threads for i in ProgressBar(1:nb_to_use)
        # identify thread and acess the MAGEMin_db of the thread
        id          = Threads.threadid()
        # println("   Thread ",id," is working on constraint ",i," of ",nb_to_use)
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        # # Calculate w_g
        w_g = calculate_w_g(variables_optim_local, variables_optim_coordinates, constraints[i].pressure, constraints[i].temperature, JOB)

        # call the forward module
        out = forward_call(JOB.solid_solution, JOB.thermodynamic_database, constraints[i], w_g, sys_in, gv, z_b, DB, splx_data)

        # check if the mineral is predicted to be stable
        if !(JOB.solid_solution in out.ph)
            # println("   Achtung: ",JOB.solid_solution," not predicted to be stable at P = $(constraints[i].pressure) kbar and T = $(constraints[i].temperature) C")
            residual_i[i] = 100
            qcmp_all[i] = 0
        else
            comp_structural_formula_clean, oxides = calc_structural_formula_element_from_output(out,JOB.solid_solution,12)

            constraint_composition = constraints[i].mineral_composition
            constraint_element = constraints[i].mineral_elements

            comp_structural_formula_clean_ordered = fix_order_structural_formula(comp_structural_formula_clean, oxides, constraint_element)

            constraint_uncertainties = bingo_generate_fake_uncertainties(constraint_composition)

            qcmp_phase = bingo_calculate_qcmp_phase(comp_structural_formula_clean_ordered,constraint_composition,constraint_uncertainties)

            # println("\n\n",qcmp_phase,"\n\n")

            residual_i[i] = (100-qcmp_phase)^2
            qcmp_all[i] = qcmp_phase
        end

    end

    # println("residual_i = ", 100 .- sqrt.(residual_i))
    # println("q_cpm      = ", qcmp_all)
    println("\n   Residual = ", sum(residual_i))
    println("   Metrics = ", sum(qcmp_all)/length(qcmp_all))

    println("\n   Margules = ", variables_optim_local)

    return sum(residual_i)
end


# """
#     objective_function(*args)

# Add an objectove function here.
# """
# function objective_function(x0, norm, JOB, constraints)
#     mineral_compositions = zeros(size(constraints)[1], 11)
#     mineral_composition = zeros(11)
#     for i = eachindex(constraints)
#         constrain = constraints[i]
#         out = forward_call(phase, database, constrain, margules, sys_in)

#         phase_idx = findfirst(x->x==phase, out.ph)
#         # extract the mineral composition in wt%
#         mineral_composition .= out.SS_vec[phase_idx].Comp_wt
#         mineral_compositions[i, :] .= mineral_composition
#     end

#     return mineral_compositions
# end

# # # test objective_function
# # database = "mp"

# # constraint_A = Constraint(8.0, 650, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])
# # constraint_B = Constraint(8.0, 600, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])

# # constraints = [constraint_A, constraint_B]

# # phase = "bi"
# # margules = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

# # @btime objective_function(margules, constraints, database, phase, "wt")
# # # 213.704 ms (896 allocations: 58.55 KiB)










# # test forward_call
# database = "mp"
# global gv, z_b, DB, splx_data  = init_MAGEMin(database);

# constraint = Constraint(8.0, 650, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])
# phase = "bi"
# margules = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

# out = forward_call(phase, constraint, gv, z_b, DB, splx_data, margules, "wt")


# """
#     objective_function(*args)

# Add an objectove function here.
# """
# function objective_function(margules, constraints, database="mp", phase = "bi", sys_in = "mol")
#     mineral_compositions = zeros(size(constraints)[1], 11)
#     mineral_composition = zeros(11)
#     for i = eachindex(constraints)
#         constrain = constraints[i]
#         out = forward_call(phase, database, constrain, margules, sys_in)

#         phase_idx = findfirst(x->x==phase, out.ph)
#         # extract the mineral composition in wt%
#         mineral_composition .= out.SS_vec[phase_idx].Comp_wt
#         mineral_compositions[i, :] .= mineral_composition
#     end

#     return mineral_compositions
# end

# # # test objective_function
# # database = "mp"

# # constraint_A = Constraint(8.0, 650, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])
# # constraint_B = Constraint(8.0, 600, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])

# # constraints = [constraint_A, constraint_B]

# # phase = "bi"
# # margules = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

# # @btime objective_function(margules, constraints, database, phase, "wt")
# # # 213.704 ms (896 allocations: 58.55 KiB)

