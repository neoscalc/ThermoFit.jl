function inversion_run(JOB, constraints)

    # get the variables to variables_optim_coordinates
    variables_optim, variables_optim_bounds, variables_optim_coordinates = get_variables_optim(JOB)

    # println(variables_optim)

    # Set parameters for the inversion
    norm = variables_optim
    x0 = variables_optim ./ norm
    nb_constraints = length(constraints)

    # # test call the objective function with the initial values (KEEP IT FOR TEST PURPOSES)
    # residual = objective_function(x0, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates)


    # perform inversion
    res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates), x0, NelderMead())

    return res, norm
end

"""
    objective_function(*args)

Add an objective function Doc here.
"""
function objective_function(x0, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, sys_in = "mol")

    variables_optim_local = x0 .* norm

    # Check if all parameters are within the bounds:
    for i = eachindex(variables_optim_local)
        if variables_optim_local[i] < variables_optim_bounds[i,1]
            return 1e20
        elseif variables_optim_local[i] > variables_optim_bounds[i,2]
            return 1e20
        end
    end

    # Initiate MAGEMin (doesn't work)
    global gv, z_b, DB, splx_data  = init_MAGEMin(JOB.thermodynamic_database);

    gv.verbose = -1

    if isequal(JOB.number_constraints_max,0)
        nb_to_use = nb_constraints
    elseif JOB.number_constraints_max > nb_constraints
        nb_to_use = nb_constraints
    else
        nb_to_use = JOB.number_constraints_max
    end

    step_print = nb_to_use ÷ 10;
    count = 0;

    println("\n-> New iteration with ",nb_to_use," constraints <-")

    residual_i = zeros(nb_to_use)
    qcmp_all = zeros(nb_to_use)
    for i = 1:nb_to_use #nb_constraints

        # print for large job
        count = count + 1
        if count == step_print
            println("      $(i)/$(nb_to_use) (",i/nb_to_use*100,"%)")
            count = 0
        end

        # Calculate w_g
        w_g = calculate_w_g(variables_optim_local, variables_optim_coordinates, constraints[i].pressure, constraints[i].temperature, JOB)

        # println(w_g)

        # call the forward module
        out = forward_call(JOB.solid_solution, JOB.thermodynamic_database, constraints[i], w_g, sys_in, gv, z_b, DB, splx_data)

        # check if the mineral is predicted to be stable
        if !(JOB.solid_solution in out.ph)
            println("   Achtung: ",JOB.solid_solution," not predicted to be stable at P = $(constraints[i].pressure) kbar and T = $(constraints[i].temperature) C")
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

    finalize_MAGEMin(gv,DB, z_b)

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

