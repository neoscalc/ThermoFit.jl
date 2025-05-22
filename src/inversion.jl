"""
    inversion.jl

This module contains functions for the inversion of thermodynamic parameters using MAGEMin.

## Structures
- `JOB`
## Functions
- `variable_optimised()`
- `inversion()`
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

    # Random seed, to be used uniformly within a single job
    # (e.g., for reproducible random starting guesses)
    random_seed               ::Union{T4, Nothing}
    
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
                 random_seed                ::Union{T4, Nothing} = nothing,
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
                                random_seed,
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
    inversion(job, constraints)
Performs the inversion of thermodynamic parameters using MAGEMin.

# Arguments
- `job::JOB`:: JOB struct containing the inversion parameters
- `constraints::Vector{Constraint}`: Vector of constraints
"""
function inversion(job, constraints; loss_f::Function=loss_Qfactor, func_rel=false, log_convergence=false)
    if log_convergence
        date = Dates.format(Dates.now(), "yyyymmdd_HHMM")
        progress_log_io = open("convergence_log_$date.log", "w")
    else
        progress_log_io = nothing
    end

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
        if !func_rel
            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f, progress_log_io=progress_log_io),
                           x0, NelderMead(),
                           Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))
        else func_rel
            res = optimize(x -> objective_function_func_relation(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f, progress_log_io=progress_log_io),
                           x0, NelderMead(),
                           Optim.Options(time_limit = max_time_seconds, iterations = max_iterations, store_trace = true))
        end
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
        res_vec = zeros(n_rand_strating_guesses)

        # loop over the random starting points and perform the optimisation using NelderMead
        for i in ProgressBar(1:n_rand_strating_guesses)
            if !func_rel
                res_i = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f, progress_log_io=progress_log_io),
                                 x0[:, i], NelderMead(),
                                 Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))

            else func_rel
                res_i = optimize(x -> objective_function_func_relation(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f, progress_log_io=progress_log_io),
                                 x0[:, i], NelderMead(),
                                 Optim.Options(time_limit = max_time_seconds, iterations = max_iterations, store_trace = true))
            end

            x_optim[:, i] = Optim.minimizer(res_i)
            res_vec[i] = Optim.minimum(res_i)
        end
        
        return x_optim, res_vec, norm

    elseif job.algorithm == "ParticleSwarm"
        # check again job.normlization as ParticleSwarm struct is initialised with/without normalised bounds
        if job.normalization == true
            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f, progress_log_io=progress_log_io),
                           x0, ParticleSwarm(; lower = optim_bounds_normed[:,1], upper = optim_bounds_normed[:,2]),
                           Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))  # n_particles = 0
        else
            res = optimize(x -> objective_function(x, job, constraints, nb_constraints, MAGEMin_db, loss_f=loss_f, progress_log_io=progress_log_io),
                           x0,
                           ParticleSwarm(; lower = optim_bounds[:,1], upper = optim_bounds[:,2]),
                           Optim.Options(time_limit = max_time_seconds, iterations = max_iterations))  # n_particles = 0
        end
        return res, norm
    else
        res = 1e20
        println("Error: algorithm not recognised")
    end

    if log_convergence
        close(progress_log_io)
    end
    Finalize_MAGEMin(MAGEMin_db)
end
