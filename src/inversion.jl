
"""
    JOB(thermodynamic_database, phase_to_be_optimised, w_names, w_initial_values, w_lower_bounds, w_upper_bounds;
        algorithm = "NelderMead", number_iterations_max = 1000, normalization = true, number_constraints_max = nothing, max_time_seconds = 300)
"""
struct JOB{T1, T2, T3, T4, T5}
    thermodynamic_database  ::T1                  # Thermodynamic database to use (mp, alk, ig, igd, mb, um)
    phase_to_be_optimised   ::T1                  # Phase to optimise (MAGEMin name)

    # THERMODYNAMIC VARIABLES TO BE OPTIMISED
    # Margules parameters: WH, WS, WV
    w_names                 ::T2
    w_initial_values        ::T3                  # (n_w x 3; columns: WH, WS, WV)
    w_lower_bounds          ::T3
    w_upper_bounds          ::T3
    # Here more thermodyn parms could be added... e.g., DQFs

    # INTERNALLY CREATED VARIABLES for the inversion
    var_optim           ::Vector{Float64}
    var_optim_bounds    ::Matrix{Float64}                       # (n x 2; columns: lower_bounds, upper_bounds)
    var_optim_names     ::Vector{String}                        # Names of the variables
    var_optim_norm      ::Vector{Float64}                       # Normalization factor
    var_optim_type      ::Vector{String}                        # Type of the variable, e.g., "WH", "WS", "WV", "DQF"... -> important for Margules!
    # For each type of variable to be optimised,
    # the coordinates are stored in a seperate vector.
    margules_optim_coord                                        # Coordinates of the Margules optimised, Vector{CartesianIndices}
    # Here more thermodyn parms could be added... e.g., DQFs

    # INVERSION PARAMETERS
    algorithm               ::T1                 
    number_iterations_max   ::T4                 
    normalization           ::T5                 
    number_constraints_max  ::Union{T4, Nothing} 
    max_time_seconds        ::T4                

    function JOB(thermodynamic_database     ::T1,
                 phase_to_be_optimised      ::T1,
                 w_names                    ::T2,
                 w_initial_values           ::T3,
                 w_lower_bounds             ::T3,
                 w_upper_bounds             ::T3;
                 algorithm                  ::T1 = "NelderMead",
                 number_iterations_max      ::T4 = 1000,
                 normalization              ::T5 = true,
                 number_constraints_max     ::Union{T4, Nothing} = nothing,
                 max_time_seconds           ::T4 = 300) where {T1<:AbstractString, T2<:AbstractArray{String},
                                                               T3<:AbstractArray{<:Real}, T4<:Integer, T5<:Bool}

        # Check if the database is existent
        if !(thermodynamic_database in ["mp", "alk", "ig", "igd", "mb", "um"])
            @error "Error: Thermodynamic database name not recognised"
        end

        # Check if the algorithm is implemented
        if !(algorithm in ["NelderMead", "ParticleSwarm"])
            @error "Error: Algorithm not recognised"
        end

        # CREATE INTERNAL VARIABLES
        # (1) Margules parameters
        margules_optim, margules_optim_bounds, margules_optim_names, margules_optim_coord, margules_n = variable_optimised(w_initial_values, w_lower_bounds, w_upper_bounds, w_names)
        margules_optim_type = repeat(["W"], margules_n)
    
        # (2) Other thermodynamic parameters (e.g., DQFs)
            ### ADD HERE
        
        var_optim           = margules_optim
        var_optim_bounds    = margules_optim_bounds    
        var_optim_names     = margules_optim_names    
        var_optim_norm      = margules_optim .+ 0.1            # use the initial values as normalization factor + "small offset"
        var_optim_type      = margules_optim_type

        new{T1, T2, T3, T4, T5}(thermodynamic_database, phase_to_be_optimised, w_names, w_initial_values, w_lower_bounds, w_upper_bounds,
                                var_optim, var_optim_bounds, var_optim_names, var_optim_norm, var_optim_type, margules_optim_coord,
                                algorithm, number_iterations_max, normalization, number_constraints_max, max_time_seconds)
    end
end


"""
    variable_optimised(lower_bounds::AbstractArray, upper_bounds::AbstractArray)
Helper func used in constructor of JOB to identify which variables are allowed to be optimised.
"""
function variable_optimised(variable::AbstractArray, lower_bounds::AbstractArray, upper_bounds::AbstractArray, variable_names::AbstractArray)
    # Variable to optimise are defined by differing lower and upper bounds
    # Check for which bounds this is the case
    mask = lower_bounds .< upper_bounds
    n = sum(mask)

    var_opti = variable[mask]
    var_opti_bounds = [lower_bounds[mask] upper_bounds[mask]]
    var_opti_names = variable_names[any(mask, dims=2)]

    var_opti_coord = findall(mask)

    return var_opti, var_opti_bounds, var_opti_names, var_opti_coord, n
end


function job_check_consistency(job)
    println("\n        *** CHECKING job CONSISTENCY ***\n")
    nb_wg = length(job.w_names)
    if isequal(size(job.w_initial_values), (nb_wg, 3)) == false
        error("w_initial_values must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_initial_values (size):        ok")
    end
    if isequal(size(job.w_lower_bounds), (nb_wg, 3)) == false
        error("w_lower_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_lower_bounds (size):          ok")
    end
    if isequal(size(job.w_upper_bounds), (nb_wg, 3)) == false
        error("w_upper_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_upper_bounds (size):          ok")
    end
    
    println("   -----------------------------------------------------")
    println(" - Variables to be optimized [name type start min max]:")
    type_w = ["WH","WS","WV"]
    for i = 1:nb_wg
        for j = 1:3
            if job.w_upper_bounds[i,j] > job.w_lower_bounds[i,j]
                println("    ", job.w_names[i], "  \t", type_w[j], "\t ",job. w_initial_values[i,j], "\t ", job.w_lower_bounds[i,j], " \t ", job.w_upper_bounds[i,j]) 
            end
        end
    end
    println("   -----------------------------------------------------")
    println(" - Algorithm:                      ", job.algorithm)
    println(" - Maximum number of iterations:   ", job.number_iterations_max)
    println(" - Normalization:                  ", job.normalization)
    println(" - Maximum number of constraints:  ", job.number_constraints_max)
end


function inversion_run(job, constraints)
    # Set parameters for the inversion
    if job.normalization == true
        x0 = job.var_optim ./ job.var_optim_norm
    else
        x0 = job.var_optim
    end
    nb_constraints = length(constraints)

    # Initiate MAGEMin
    MAGEMin_db  = Initialize_MAGEMin(job.thermodynamic_database, verbose=false);

    # # test call the objective function with the initial values (KEEP IT FOR TEST PURPOSES)
    # residual = objective_function(x0, job, constraints, nb_constraints)
    max_time_seconds = job.max_time_seconds

    # perform inversion using the Nelder-Mead algorithm
    if job.algorithm == "NelderMead"
        res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db), x0, NelderMead(), Optim.Options(time_limit = max_time_seconds, iterations = job.number_iterations_max))
    elseif job.algorithm == "ParticleSwarm"

        if job.normalization == true
            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db), x0, ParticleSwarm(; lower = job.var_optim_bounds[:,1] ./ jon.var_optim_norm, upper = job.var_optim_bounds[:,2] ./ jon.var_optim_norm), Optim.Options(time_limit = max_time_seconds, iterations = job.number_iterations_max))  # n_particles = 0
        else
            println("x0 = ", x0)
            println("lower = ", [job.var_optim_bounds[:,1]])
            println("upper = ", [job.var_optim_bounds[:,2]])
            # overwrite norm to be 1
            norm = ones(length(x0))

            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db), x0, ParticleSwarm(; lower = job.var_optim_bounds[:,1], upper = job.var_optim_bounds[:,2]), Optim.Options(time_limit = max_time_seconds, iterations = job.number_iterations_max))  # n_particles = 0
        end
    else
        res = 1e20
        println("Error: algorithm not recognised")
    end

    Finalize_MAGEMin(MAGEMin_db)

    return res, job.var_optim_norm
end

"""
    objective_function(*args)

Add an objective function Doc here.
"""
function objective_function(x0, job, constraints, nb_constraints, MAGEMin_db, sys_in = "mol")

    # define number of constraints to use for the inversion
    if isequal(job.number_constraints_max,0)
        nb_to_use = nb_constraints
    elseif job.number_constraints_max > nb_constraints
        nb_to_use = nb_constraints
    else
        nb_to_use = job.number_constraints_max
    end

    println("\n-> New iteration with ",nb_to_use," constraints <-")

    # Denormalise variables to optimise (Margules) for G-minimisation
    if job.normalization == true
        variables_optim_local = x0 .* job.var_optim_norm
    else
        variables_optim_local = x0
    end

    # Check if all parameters are within the bounds, if not return a very high residual:
    for i = eachindex(variables_optim_local)
        if variables_optim_local[i] < job.var_optim_bounds[i,1]
            return 1e20
        elseif variables_optim_local[i] > job.var_optim_bounds[i,2]
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

        # Calculate w_g
        w_optimised     = variables_optim_local[job.var_optim_type .== "W"]
        P_kbar          = constraints[i].pressure_GPa * 10
        T_K             = constraints[i].temperature_C + 273.15
        w_matrix        = copy(job.w_initial_values)   # copy needed for multi_threading, so every thread has its own copy of w_all
        w_coordinates   = job.margules_optim_coord
        
        w_g = calculate_w_g(w_optimised, P_kbar, T_K, w_matrix, w_coordinates)

        # call the forward module
        out = forward_call(job.phase_to_be_optimised, job.thermodynamic_database, constraints[i], w_g, sys_in, gv, z_b, DB, splx_data)

        # check if the mineral is predicted to be stable
        if !(job.phase_to_be_optimised in out.ph)
            # println("   Achtung: ",job.phase_to_be_optimised," not predicted to be stable at P = $(constraints[i].pressure) kbar and T = $(constraints[i].temperature) C")
            residual_i[i] = 100
            qcmp_all[i] = 0
        else
            composition_predicted = out.SS_vec[findfirst(x->x==job.phase_to_be_optimised, out.ph)].Comp_apfu
            #reorder to match the order of elements in the constraint
            idx_elements_constraint_in_out = indexin(constraints[i].mineral_elements, out.elements)
            if nothing in idx_elements_constraint_in_out
                @error "Error: Constraint #$(i) has a mineral composition with elements that are not in the predicted mineral composition."
            else
                composition_predicted = composition_predicted[idx_elements_constraint_in_out]
            end

            # compare predicted composition with the constraint composition
            # calculate a loss (q_cmp)
            constraint_composition = constraints[i].mineral_composition_apfu[job.phase_to_be_optimised]
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
