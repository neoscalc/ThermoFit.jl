"""
    inversion.jl

This module contains functions for the inversion of thermodynamic parameters using MAGEMin.

## Structures
- `JOB`
## Functions
- `variable_optimised()`
- `print_job()`
- `inversion()`
- `objective_function()`
"""


"""
    JOB struct to set up the inversion:
- THERMODYNAMIC VARIABLES TO BE OPTIMISED
- INVERSION PARAMETERS
"""
struct JOB{T1, T2, T3, T4, T5, T6}
    thermodynamic_database  ::T1                  # Thermodynamic database to use (mp, alk, ig, igd, mb, um)
    phase_to_be_optimised   ::T1                  # Phase to optimise (MAGEMin name)

    # THERMODYNAMIC VARIABLES TO BE OPTIMISED
    # Margules parameters: WH, WS, WV
    w_names                 ::Union{T2, Nothing}
    w_initial_values        ::Union{T3, Nothing}                  # (n_w x 3; columns: WH, WS, WV)
    w_lower_bounds          ::Union{T3, Nothing}
    w_upper_bounds          ::Union{T3, Nothing}
    # G correction for endmembers ("DQFs")
    g0_corr_endmembers      ::Union{T2, Nothing}
    g0_corr_initial_values  ::Union{T6, Nothing}
    g0_corr_lower_bounds    ::Union{T6, Nothing}
    g0_corr_upper_bounds    ::Union{T6, Nothing}
    # Here more thermodyn parms could be added...

    # INTERNALLY CREATED VARIABLES for the inversion
    var_optim           ::Vector{Float64}
    var_optim_bounds    ::Matrix{Float64}                       # (n x 2; columns: lower_bounds, upper_bounds)
    var_optim_names     ::Vector{String}                        # Names of the variables
    var_optim_norm      ::Vector{Float64}                       # Normalization factor
    var_optim_type      ::Vector{String}                        # Type of the variable, e.g., "WH", "WS", "WV", "G0_corr"... -> important for Margules!
    # For each type of variable to be optimised,
    # the coordinates are stored in a seperate vector.
    margules_optim_coord ::Union{Vector{CartesianIndex{2}}, Nothing}       # Coordinates of the Margules optimised, Vector{CartesianIndices}
    g0_corr_optim_coord  ::Union{Vector{Int64}, Nothing}                   # Coordinates of the G-correction optimised, Vector{CartesianIndices}
    # Here more thermodyn parms could be added...
    
    # INVERSION PARAMETERS
    algorithm                 ::T1                 
    number_iterations_max     ::T4                 
    normalization             ::T5                 
    number_constraints_max    ::Union{T4, Nothing} 
    max_time_seconds          ::T4
    n_rand_strating_guesses   ::Union{T4, Nothing}
    
    # VERBOSITY
    verbose                   ::T5

    function JOB(thermodynamic_database     ::T1,
                 phase_to_be_optimised      ::T1;
                 w_names                    ::Union{T2, Nothing} = nothing,
                 w_initial_values           ::Union{T3, Nothing} = nothing,
                 w_lower_bounds             ::Union{T3, Nothing} = nothing,
                 w_upper_bounds             ::Union{T3, Nothing} = nothing,
                 g0_corr_endmembers         ::Union{T2, Nothing} = nothing,
                 g0_corr_initial_values     ::Union{T6, Nothing} = nothing,
                 g0_corr_lower_bounds       ::Union{T6, Nothing} = nothing,
                 g0_corr_upper_bounds       ::Union{T6, Nothing} = nothing,
                 algorithm                  ::T1                 = "NelderMead",
                 number_iterations_max      ::T4                 = 1000,
                 normalization              ::T5                 = true,
                 number_constraints_max     ::Union{T4, Nothing} = nothing,
                 max_time_seconds           ::T4                 = 300,
                 n_rand_strating_guesses    ::Union{T4, Nothing} = nothing,
                 verbose                    ::T5                 = true) where {T1<:AbstractString, T2<:AbstractArray{String},
                                                                                T3, T4<:Integer, T5<:Bool, T6}

        # Check if the database is existent
        if !(thermodynamic_database in ["mp", "alk", "ig", "igd", "mb", "um"])
            @error "Error: Thermodynamic database name not recognised"
        end

        # Check if the algorithm is implemented
        if !(algorithm in ["NelderMead", "NelderMead_random_guess", "ParticleSwarm"])
            @error "Error: Algorithm not recognised"
        end

        # Check if the thermodynmaic parameters to be optimised are consistent
        # (1) Margules parameters
        if !isnothing(w_names)
            nb_wg = length(w_names)
            if isequal(size(w_initial_values), (nb_wg, 3)) == false
                @error("Error: w_initial_values must be a matrix of size (nb_wg, 3)")
            end
            if isequal(size(w_lower_bounds), (nb_wg, 3)) == false
                @error("Error: w_lower_bounds must be a matrix of size (nb_wg, 3)")
            end
            if isequal(size(w_upper_bounds), (nb_wg, 3)) == false
                @error("Error: w_upper_bounds must be a matrix of size (nb_wg, 3)")
            end
        end
        # (2) G-correction parameters
        if !isnothing(g0_corr_endmembers)
            nb_g0 = length(g0_corr_endmembers)
            if isequal(size(g0_corr_initial_values), (nb_g0, )) == false
                @error("Error: g0_corr_initial_values must be a vector of size (nb_g0,)")
            end
        end

        # CREATE INTERNAL VARIABLES
        # (1) Margules parameters
        if !isnothing(w_names)
            margules_optim, margules_optim_bounds, margules_optim_names, margules_optim_coord, n_margules = variable_optimised(w_initial_values, w_lower_bounds, w_upper_bounds, w_names)
            margules_optim_type = repeat(["W"], n_margules)
        else
            n_margules = 0
            margules_optim_coord = nothing
        end
    
        # (2) G-correction parameters
        if !isnothing(g0_corr_endmembers)
            g0_corr_optim, g0_corr_optim_bounds, g0_corr_optim_names, g0_corr_optim_coord, n_g0_corr = variable_optimised(g0_corr_initial_values, g0_corr_lower_bounds, g0_corr_upper_bounds, g0_corr_endmembers)
            g0_optim_type = repeat(["G0_corr"], n_g0_corr)
        else
            n_g0_corr = 0
            g0_corr_optim_coord = nothing
        end

        if n_margules == 0 && n_g0_corr == 0
            @error "Error: No thermodynamic parameters to be optimised"
        elseif n_margules == 0
            var_optim           = g0_corr_optim
            var_optim_bounds    = g0_corr_optim_bounds
            var_optim_names     = g0_corr_optim_names
            var_optim_norm      = abs.(g0_corr_optim) .+ 0.01       # use the initial values as normalization factor + "small offset"
            var_optim_type      = g0_optim_type
        elseif n_g0_corr == 0
            var_optim           = margules_optim
            var_optim_bounds    = margules_optim_bounds
            var_optim_names     = margules_optim_names
            var_optim_norm      = abs.(margules_optim) .+ 0.01       # use the initial values as normalization factor + "small offset"
            var_optim_type      = margules_optim_type
        else
            var_optim           = vcat(margules_optim, g0_corr_optim)
            var_optim_bounds    = vcat(margules_optim_bounds, g0_corr_optim_bounds)
            var_optim_names     = vcat(margules_optim_names, g0_corr_optim_names)
            var_optim_norm      = vcat(abs.(margules_optim), abs.(g0_corr_optim)) .+ 0.01       # use the initial values as normalization factor + "small offset"
            var_optim_type      = vcat(margules_optim_type, g0_optim_type)
        end
        
        #//FIXME -  - Typing might needs a cleaner solution here, instead of dirty fix
        U1 = typeof(w_initial_values)
        U2 = typeof(g0_corr_initial_values)

        new{T1, T2, U1, T4, T5, U2}(thermodynamic_database,
                                phase_to_be_optimised,
                                w_names, w_initial_values, w_lower_bounds, w_upper_bounds,
                                g0_corr_endmembers, g0_corr_initial_values, g0_corr_lower_bounds, g0_corr_upper_bounds,
                                var_optim,
                                var_optim_bounds,
                                var_optim_names,
                                var_optim_norm,
                                var_optim_type,
                                margules_optim_coord,
                                g0_corr_optim_coord,
                                algorithm,
                                number_iterations_max,
                                normalization,
                                number_constraints_max,
                                max_time_seconds,
                                n_rand_strating_guesses,
                                verbose)
    end
end


"""
    variable_optimised(lower_bounds::AbstractArray, upper_bounds::AbstractArray)
Auxillary func used in constructor of JOB to identify which variables are allowed to be optimised.
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


"""
    print_job(job)
Prints the job parameters.
"""
function print_job(job)
    println("   -----------------------------------------------------")
    println(" - Margules to be optimized [name type start min max]:")
    type_w = ["WH","WS","WV"]
    for i = 1:length(job.w_names)
        for j = 1:3
            if job.w_upper_bounds[i,j] > job.w_lower_bounds[i,j]
                println("    ", job.w_names[i],
                        "  \t", type_w[j],
                        "\t ",job. w_initial_values[i,j],
                        "\t ", job.w_lower_bounds[i,j],
                        " \t ", job.w_upper_bounds[i,j]) 
            end
        end
    end
    if !isnothing(job.g0_corr_endmembers)
        println(" - G0 corrections [name start min max]:")
        for i = 1:length(job.g0_corr_endmembers)
            if job.g0_corr_upper_bounds[i] > job.g0_corr_lower_bounds[i]
                println("    ", job.g0_corr_endmembers[i],
                        "  \t", job.g0_corr_initial_values[i],
                        "\t ", job.g0_corr_lower_bounds[i],
                        " \t ", job.g0_corr_upper_bounds[i])
            end
        end
    end
    println("   -----------------------------------------------------")
    println(" - Algorithm:                      ", job.algorithm)
    println(" - Maximum number of iterations:   ", job.number_iterations_max)
    println(" - Normalization:                  ", job.normalization)
    println(" - Maximum number of constraints:  ", job.number_constraints_max)
end


"""
    inversion(job, constraints)
Performs the inversion of thermodynamic parameters using MAGEMin.

# Arguments
- `job::JOB`:: JOB struct containing the inversion parameters
- `constraints::Vector{Constraint}`: Vector of constraints
"""
function inversion(job, constraints; loss_f::Function=loss_Qfactor)
    # Unpack variables form JOB struct here
    var_optim               = job.var_optim
    norm                    = job.var_optim_norm
    optim_bounds            = job.var_optim_bounds
    optim_bounds_normed     = optim_bounds ./ norm
    
    if job.normalization == true
        x0 = var_optim ./ norm
    else
        x0 = var_optim
        norm = ones(length(x0))
    end

    algorithm = job.algorithm
    max_time_seconds = job.max_time_seconds
    max_iterations = job.number_iterations_max
    n_rand_strating_guesses = job.n_rand_strating_guesses

    # determine the number of constraints to use for the inversion
    nb_constraints = length(constraints)
    # check if the number of constraints is limited by the user in JOB
    if !isnothing(job.number_constraints_max)
        if job.number_constraints_max < nb_constraints
            nb_constraints = job.number_constraints_max
        end
    end
    if job.verbose
        println("   Number of constraints used for the inversion: ", nb_constraints)
    end


    # Initiate MAGEMin
    MAGEMin_db  = Initialize_MAGEMin(job.thermodynamic_database, solver=2, verbose=false);


    # PERFORM INVERSION: Check which algorithm to use
    if algorithm == "NelderMead"
        res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f),
                       x0, NelderMead(),
                       Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))
        return res, norm

    elseif algorithm == "NelderMead_random_guess"
        # create a vector of random starting guesses for the 'var_optim' within the 'var_opti_bounds'
        var_optim_starting_points = rand(length(var_optim), n_rand_strating_guesses) .* (optim_bounds[:,2] .- optim_bounds[:,1]) .+ optim_bounds[:,1]
        
        if job.normalization == true
            x0 = var_optim_starting_points ./ norm
        else
            x0 = var_optim_starting_points
            norm = ones(size(x0, 1))
        end

        # create a vector for res of the optimisation
        x_optim = zeros(size(x0))
        residual_vec = zeros(n_rand_strating_guesses)

        # loop over the random starting points and perform the optimisation using NelderMead
        for i in ProgressBar(1:n_rand_strating_guesses)
            res_i = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f),
                             x0[:, i], NelderMead(),
                             Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))
        
            x_optim[:, i] = Optim.minimizer(res_i)
            residual_vec[i] = Optim.minimum(res_i)
        end
        
        return x_optim, residual_vec, norm

    elseif job.algorithm == "ParticleSwarm"
        # check again job.normlization as ParticleSwarm struct is initialised with/without normalised bounds
        if job.normalization == true
            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f),
                           x0, ParticleSwarm(; lower = optim_bounds_normed[:,1], upper = optim_bounds_normed[:,2]),
                           Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))  # n_particles = 0
        else
            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f),
                           x0,
                           ParticleSwarm(; lower = optim_bounds[:,1], upper = optim_bounds[:,2]),
                           Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))  # n_particles = 0
        end
        return res, norm
    else
        res = 1e20
        println("Error: algorithm not recognised")
    end

    Finalize_MAGEMin(MAGEMin_db)
end


"""
    objective_function(*args)

Add an objective function Doc here.
"""
function objective_function(x0, job, constraints, nb_constraints, MAGEMin_db; loss_f::Function=loss_Qfactor)
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

    # Initiate vectors for loss (to be minimised) and Q_cmp (used as metric that is printed during the inversion)
    loss_vec            = zeros(nb_constraints)
    qcmp_vec            = zeros(nb_constraints)
    phase_pred_stable   = zeros(nb_constraints)

    # Conditional definition of the iterator before the looping over the constraints,
    # this way the verbosity can be toggled, by defining the iterator either with or
    # without a progressbar
    if job.verbose
        println("\n-> New iteration <-")
        itr = ProgressBar(1:nb_constraints)
    else
        itr = 1:nb_constraints
    end

    # --------------------------------------------------
    # MAIN LOOP: Calculate the loss over all constraints
    # --------------------------------------------------
    @threads for i in itr
        constraint = constraints[i]

        constraint, gv, z_b, DB, splx_data, w_g, g0_corr, g0_corr_endmembers = prepare_forward_call(variables_optim_local,
                                                                                                    MAGEMin_db,
                                                                                                    job, 
                                                                                                    constraint)

        # call the forward module
        out = forward_call(job.phase_to_be_optimised,
                           job.thermodynamic_database,
                           constraint,
                           gv, z_b, DB, splx_data,
                           w_g = w_g,
                           g0_corr_vec = g0_corr,
                           g0_corr_em = g0_corr_endmembers)

        # --------------------------------
        # Process the minimisation output:
        # --------------------------------
        # check if the mineral is predicted to be stable
        if !(job.phase_to_be_optimised in out.ph)
            # println("   Achtung: ",job.phase_to_be_optimised," not predicted to be stable at P = $(constraints[i].pressure_GPa) kbar and T = $(constraints[i].temperature_C) C")
            qcmp_vec[i] = 0
            loss_vec[i] = 100
        else
            # change 0 > 1 in the phase_pred_stable vector
            phase_pred_stable[i] = 1

            composition_predicted = out.SS_vec[findfirst(x->x==job.phase_to_be_optimised, out.ph)].Comp_apfu
            #reorder the compositional vector to match the order of elements in the constraint
            idx_elements_constraint_in_out = indexin(constraints[i].mineral_elements, out.elements)
            if nothing in idx_elements_constraint_in_out
                @error "Error: Constraint #$(i) has a mineral composition with elements that are not in the predicted mineral composition."
            else
                composition_predicted = composition_predicted[idx_elements_constraint_in_out]
            end

            constraint_composition = constraints[i].mineral_composition_apfu[job.phase_to_be_optimised]
            # calculate loss: Default is loss_Qfactor: (100 - Q_cmp)
            loss_vec[i] = loss_f(composition_predicted, constraint_composition)
            # calculate Q_cmp as a metric beside the loss function
            qcmp_vec[i] = quality_factor(constraint_composition, composition_predicted)            
        end

    end

    # calculate the mean loss over all constraints and the fraction of constraints where the phase optimised is predicted stable
    mean_loss = sum(loss_vec) / nb_constraints
    frac_phase_present = (sum(phase_pred_stable) / nb_constraints) * 100
    # calculate the mean between the two
    residual = .5 * (mean_loss + frac_phase_present)

    # Calculate the mean Q_cmp metric
    qcmp = sum(qcmp_vec) / nb_constraints

    if job.verbose
        println("\n   Residual = ", residual)
        println("   Metrics = ", qcmp)
        println("   Fraction of constraints where the phase is predicted stable = ", frac_phase_present)
        println("   Loss composition = ", mean_loss)
        println("\n   Optimied variables = ", variables_optim_local)
    end

    return residual
end


"""
    objective_function(*args)

Add an objective function Doc here.
"""
function objective_function_func_relation(x0, job, constraints, nb_constraints, MAGEMin_db; loss_f::Function=Ti_sat_misfit)
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

    # Initiate vectors for the loss
    loss_vec            = zeros(nb_constraints)
    metric_vec          = zeros(nb_constraints)

    # Conditional definition of the iterator before the looping over the constraints,
    # this way the verbosity can be toggled, by defining the iterator either with or
    # without a progressbar
    if job.verbose
        println("\n-> New iteration <-")
        itr = ProgressBar(1:nb_constraints)
    else
        itr = 1:nb_constraints
    end

    # --------------------------------------------------
    # MAIN LOOP: Calculate the loss over all constraints
    # --------------------------------------------------
    @threads for i in itr
        constraint = constraints[i]

        constraint, gv, z_b, DB, splx_data, w_g, g0_corr, g0_corr_endmembers = prepare_forward_call(variables_optim_local,
                                                                                                    MAGEMin_db,
                                                                                                    job, 
                                                                                                    constraint)

        # call the forward module
        out = forward_call(job.phase_to_be_optimised,
                           job.thermodynamic_database,
                           constraint,
                           gv, z_b, DB, splx_data,
                           w_g = w_g,
                           g0_corr_vec = g0_corr,
                           g0_corr_em = g0_corr_endmembers)

        # --------------------------------
        # Process the minimisation output:
        # --------------------------------
        # check if the mineral is predicted to be stable
        if !(job.phase_to_be_optimised in out.ph)
            # As the stabilty field of biotite should not be altered significantly by optimising
            # any relation defined for its composition, "losing" biotite form the assemblage
            # should be penalised by a large loss.
            # A loss equivalent to a misfit between T_model and T_Henry of 400°C is taken as appropriate.
            loss_vec[i] = 1e20
        else
            composition_predicted = out.SS_vec[findfirst(x->x==job.phase_to_be_optimised, out.ph)].Comp_apfu

            loss_vec[i] = loss_f(composition_predicted, constraint.temperature_C, out.ph, out.elements)
            metric_vec[i] = Ti_in_Bt_misfit(composition_predicted, constraint.temperature_C, out.elements)          
        end

    end

    # calculate the mean loss over all constraints and the fraction of constraints where the phase optimised is predicted stable
    residual = sum(loss_vec) / nb_constraints
    metric   = sum(skipmissing(metric_vec)) / nb_constraints

    if job.verbose
        println("\n   Residual = ", residual)
        println("   ΔT(model-Henry05) = ", metric)
        println("\n   Optimied variables = ", variables_optim_local)
    end

    return residual
end
