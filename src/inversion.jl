
"""
    JOB(thermodynamic_database, phase_to_be_optimised, w_names, w_initial_values, w_lower_bounds, w_upper_bounds;
        algorithm = "NelderMead", number_iterations_max = 1000, normalization = true, number_constraints_max = nothing, max_time_seconds = 300)
"""
struct JOB{T1, T2, T3, T4, T5}
    thermodynamic_database  ::T1                  # Thermodynamic database to use (mp, alk, ig, igd, mb, um)
    phase_to_be_optimised   ::T1                  # Phase to optimise (MAGEMin name)

    w_names                 ::T2
    w_initial_values        ::T3                  # (n_w x 3; columns: WH, WS, WV)
    w_lower_bounds          ::T3
    w_upper_bounds          ::T3

    algorithm         ::T1                 
    number_iterations_max   ::T4                 
    normalization           ::T5                 
    number_constraints_max  ::Union{T4, Nothing} 
    max_time_seconds        ::T4                

    function JOB(thermodynamic_database::T1, phase_to_be_optimised::T1, w_names::T2,
                 w_initial_values::T3, w_lower_bounds::T3, w_upper_bounds::T3;
                 algorithm::T1 = "NelderMead", number_iterations_max::T4 = 1000,
                 normalization::T5 = true, number_constraints_max ::Union{T4, Nothing} = nothing,
                 max_time_seconds::T4 = 300) where {T1<:AbstractString, T2<:AbstractArray{String},
                                                    T3<:AbstractArray{<:Real}, T4<:Integer, T5<:Bool}

        # Check if the database is existent
        if !(thermodynamic_database in ["mp", "alk", "ig", "igd", "mb", "um"])
            @error "Error: Thermodynamic database name not recognised"
        end

        # Check if the algorithm is implemented
        if !(algorithm in ["NelderMead", "ParticleSwarm"])
            @error "Error: Algorithm not recognised"
        end
        
        new{T1, T2, T3, T4, T5}(thermodynamic_database, phase_to_be_optimised, w_names, w_initial_values, w_lower_bounds, w_upper_bounds,
            algorithm, number_iterations_max, normalization, number_constraints_max, max_time_seconds)
    end
end


function inversion_run(JOB, constraints)
    # get the variables to variables_optim_coordinates
    variables_optim, variables_optim_bounds, variables_optim_coordinates = get_variables_optim(JOB)

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
        res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db), x0, NelderMead(), Optim.Options(time_limit = max_time_seconds, iterations = JOB.number_iterations_max))
    elseif JOB.algorithm == "ParticleSwarm"

        if JOB.normalization == true
            res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db), x0, ParticleSwarm(; lower = variables_optim_bounds[:,1] ./ norm, upper = variables_optim_bounds[:,2] ./ norm), Optim.Options(time_limit = max_time_seconds, iterations = JOB.number_iterations_max))  # n_particles = 0
        else
            println("x0 = ", x0)
            println("lower = ", [variables_optim_bounds[:,1]])
            println("upper = ", [variables_optim_bounds[:,2]])
            # overwrite norm to be 1
            norm = ones(length(x0))

            res = optimize(x -> objective_function(x, norm, JOB, constraints, nb_constraints, variables_optim_bounds, variables_optim_coordinates, MAGEMin_db), x0, ParticleSwarm(; lower = variables_optim_bounds[:,1], upper = variables_optim_bounds[:,2]), Optim.Options(time_limit = max_time_seconds, iterations = JOB.number_iterations_max))  # n_particles = 0
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


    @threads for i in ProgressBar(1:nb_to_use)
        # identify thread and acess the MAGEMin_db of the thread
        id          = Threads.threadid()
        # println("   Thread ",id," is working on constraint ",i," of ",nb_to_use)
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        # # Calculate w_g
        w_g = calculate_w_g(variables_optim_local, variables_optim_coordinates, constraints[i].pressure_GPa, constraints[i].temperature_C, JOB)
        # call the forward module
        out = forward_call(JOB.phase_to_be_optimised, JOB.thermodynamic_database, constraints[i], w_g, sys_in, gv, z_b, DB, splx_data)

        # check if the mineral is predicted to be stable
        if !(JOB.phase_to_be_optimised in out.ph)
            # println("   Achtung: ",JOB.phase_to_be_optimised," not predicted to be stable at P = $(constraints[i].pressure) kbar and T = $(constraints[i].temperature) C")
            residual_i[i] = 100
            qcmp_all[i] = 0
        else
            composition_predicted = out.SS_vec[findfirst(x->x==JOB.phase_to_be_optimised, out.ph)].Comp_apfu
            #reorder to match the order of elements in the constraint
            idx_elements_constraint_in_out = indexin(constraints[i].mineral_elements, out.elements)
            if nothing in idx_elements_constraint_in_out
                @error "Error: Constraint #$(i) has a mineral composition with elements that are not in the predicted mineral composition."
            else
                composition_predicted = composition_predicted[idx_elements_constraint_in_out]
            end

            # compare predicted composition with the constraint composition
            # calculate a loss (q_cmp)
            constraint_composition = constraints[i].mineral_composition_apfu[JOB.phase_to_be_optimised]
            constraint_uncertainties = bingo_generate_fake_uncertainties(constraint_composition)
            qcmp_phase = bingo_calculate_qcmp_phase(composition_predicted,constraint_composition,constraint_uncertainties)

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
