using Parameters

@with_kw struct constraint
    pressure::Float64
    temperature::Float64
    bulk::Vector{Float64}
    bulk_oxides::Vector{String}
    mineral_composition::Vector{Float64}
    mineral_elements::Vector{String}
end

function inversion_run(JOB, constraints)

    # get the variables to variables_optim_coordinates
    variables_optim, variables_optim_bounds, variables_optim_coordinates = get_variables_optim(JOB)
    
    # println(variables_optim)

    # test call the objective function with the initial values
    norm = variables_optim
    x0 = variables_optim ./ norm
    nb_constraints = length(constraints)

    residual = objective_function(x0, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates)

    println(residual)

    # perform inversion



end

"""
    objective_function(*args)

Add an objectove function here.
"""
function objective_function(x0, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates)
    
    variables_optim_local = x0 .* norm

    # Implement the check for bounds here (NOT DONE YET)!!!!!!!!!!!!

    # Initiate MAGEMin (doesn't work)
    # global gv, z_b, DB, splx_data  = init_MAGEMin(JOB.thermodynamic_database);

    residual_i = zeros(nb_constraints)

    for i = 1:1#nb_constraints
        
        # Calculate w_g
        w_g = calculate_w_g(variables_optim_local, variables_optim_coordinates, constraints[i].pressure, constraints[i].temperature, JOB)

        # println(w_g)

        # call the forward module
        out = forward_call(JOB.solid_solution, JOB.thermodynamic_database, constraints[i], w_g, "wt")

        comp_structural_formula_clean, oxides = calc_structural_formula_element_from_output(out,"bi",12)

        constraint_composition = constraints[i].mineral_composition
        constraint_element = constraints[i].mineral_elements

        comp_structural_formula_clean_ordered = fix_order_structural_formula(comp_structural_formula_clean, oxides, constraint_element)

        constraint_uncertainties = bingo_generate_fake_uncertainties(constraint_composition)

        qcmp_phase = bingo_calculate_qcmp_phase(comp_structural_formula_clean_ordered,constraint_composition,constraint_uncertainties)

        println("\n\n",qcmp_phase,"\n\n")

        residual_i[i] = (100-qcmp_phase)^2
    end

    return sum(residual_i)
end

function forward_call(phase, database, constraint, w_g, sys_in)
    global gv, z_b, DB, splx_data  = init_MAGEMin(database);

    pressure = constraint.pressure
    temperature = constraint.temperature
    bulk = constraint.bulk
    bulk_oxides = constraint.bulk_oxides

    gv = define_bulk_rock(gv, bulk, bulk_oxides, sys_in, database);
    gv, z_b, DB, splx_data = pwm_init(pressure, temperature, gv, z_b, DB, splx_data);

    ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
    ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    ss_idx = findfirst(x->x==phase, ss_names);

    W = unsafe_wrap(Vector{Cdouble},ss_struct[ss_idx].W, ss_struct[ss_idx].n_w)
    W = w_g

    out = pwm_run(gv, z_b, DB, splx_data)

    finalize_MAGEMin(gv,DB, z_b)

    return out
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

