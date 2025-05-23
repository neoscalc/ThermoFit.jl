
"""
    print_job(job)
Prints the job parameters.
"""
function print_job(job  ::JOB;
                   io   ::IO = stdout)
    date = Dates.format(Dates.now(), "dd.mm.yyyy HH:MM:SS")

    println(io, "-----------------------------------------------------")
    println(io, "              THERMOFIT INVERSION RUN")
    println(io, "-----------------------------------------------------")
    println(io, " Run started on: ", date)
    println(io, "-----------------------------------------------------")
    println(io, " INVERSION RUN SETUP:")
    println(io, " - Thermodynamic database:            ", job.thermodynamic_database)
    println(io, " - Phase to be optimised:             ", job.phase_to_be_optimised)
    println(io, " - Number of variables to be optimised:", length(job.var_optim))
    println(io, "\n")
    println(io, " - Variables to be optimised [name (type) start min max]:")
    println(io, "-----------------------------------------------------")
    println(io, "    - Margules parameters:")
    type_w = ["WH","WS","WV"]
    for i = 1:length(job.w_names)
        for j = 1:3
            if job.w_upper_bounds[i,j] > job.w_lower_bounds[i,j]
                println(io, "    ", job.w_names[i],
                        "    ", type_w[j],
                        "    ",job. w_initial_values[i,j],
                        "    ", job.w_lower_bounds[i,j],
                        "    ", job.w_upper_bounds[i,j]) 
            end
        end
    end
    if !isnothing(job.g0_corr_endmembers)
        println(io, "    - G° corrections [name start min max]:")
        for i = 1:length(job.g0_corr_endmembers)
            if job.g0_corr_upper_bounds[i] > job.g0_corr_lower_bounds[i]
                println(io, "    ", job.g0_corr_endmembers[i],
                        "    ", job.g0_corr_initial_values[i],
                        "    ", job.g0_corr_lower_bounds[i],
                        "    ", job.g0_corr_upper_bounds[i])
            end
        end
    end
    println(io, "\n")
    println(io, "Variables to be optimised normalised?: ", job.normalization)
    println(io, " - Normalisation factors: ", job.var_optim_norm)
    println(io, "\n")
    println(io, "-----------------------------------------------------")
    println(io, " - Algorithm:                      ", job.algorithm)
    println(io, " - Maximum number of iterations:   ", job.number_iterations_max)
    println(io, " - Normalization:                  ", job.normalization)
    println(io, " - Maximum number of constraints:  ", job.number_constraints_max)
    println(io, " - Maximum time (s):               ", job.max_time_seconds)
    println(io, " - Number of random starting guesses: ", job.n_rand_strating_guesses)
    println(io, " - Seed for random number generator: ", job.random_seed)
    println(io, "\n")
end


function print_constraints(nb_constraints   ::Number;
                           constraints_yaml :: Union{Nothing, String},
                           constraints_gen  ::Bool= false,
                           P_MIN_GPa        ::Union{Nothing, AbstractFloat}=nothing,
                           P_MAX_GPa        ::Union{Nothing, AbstractFloat}=nothing,
                           T_MIN_C          ::Union{Nothing, AbstractFloat}=nothing,
                           T_MAX_C          ::Union{Nothing, AbstractFloat}=nothing,
                           bulk_rocks       ::Union{Nothing, String}=nothing,
                           sys_in           ::Union{Nothing, String}=nothing,
                           io               ::IO = stdout)
    println(io, "-----------------------------------------------------")
    println(io, "          CONSTRAINTS FOR THE INVERSION")
    println(io, "-----------------------------------------------------")
    println(io, " - Number of constraints available: ", nb_constraints)
    println(io, " - Constraints file:                ", constraints_yaml)
    println(io, " - Constraints generated:           ", constraints_gen)
    println(io, "   using the following parameters:")
    println(io, "       - P_min [GPa]:           ", P_MIN_GPa)
    println(io, "       - P_max [GPa]:           ", P_MAX_GPa)
    println(io, "       - T_min [°C]:            ", T_MIN_C)
    println(io, "       - T_max [°C]:            ", T_MAX_C)
    println(io, "       - Bulk rocks:            ", bulk_rocks)
    println(io, "       - System in:             ", sys_in)
    println(io, "\n")
    println(io, "-----------------------------------------------------")
end


function print_results(optim_result     ::Optim.MultivariateOptimizationResults,
                       job              ::JOB,
                       loss_f           ::Function,
                       metric_f         ::Function;
                       io               ::IO = stdout)
    date = Dates.format(Dates.now(), "dd.mm.yyyy HH:MM:SS")

    # check what terminated the inv run and construct "conv_status" string
    converged = Optim.converged(optim_result)
    conv_status = converged ? "Converged" : "Stopped by "
    if Optim.iteration_limit_reached(optim_result)
        conv_status *= "maximum number of $(Optim.iterations(optim_result)) iterations"
    end
        if Optim.f_increased(optim_result) && !(Optim.iteration_limit_reached(optim_result))
        conv_status *= " objective increasing between iterations"
    end
    if Optim.time_run(optim_result) > Optim.time_limit(optim_result)
        conv_status *= " time limit of $(Optim.time_limit(optim_result)) seconds"
    end

    # Compute the initial and final objective function values and metrics
    # NOTE - can't be bothered atm to calculate as the objective function and metrics
    # differ between the inversion and the inversion to a functional relation
    # these data are already logged in the convergence log and can be found there
    PLACEHOLDER = "\n These values can be found in the convergence log."

    # Extract the optimised thermodynamic parameters
    opt_params = Optim.minimizer(optim_result) .* job.var_optim_norm
    opt_params_names = job.var_optim_names



    println(io, "-----------------------------------------------------")
    println(io, "             RESULTS OF THE INVERSION")
    println(io, "-----------------------------------------------------")
    println(io, " Run terminated on: ", date)
    println(io, " Run terminated by: ", conv_status)
    println(io, "-----------------------------------------------------")
    println(io, " INVERSION RUN RESULTS:")
    println(io, "- Algorithm:                      ", summary(optim_result))
    println(io, "- Loss function:                  ", string(loss_f))
    println(io, "- Metric function:                ", string(metric_f))
    println(io, "- Number of iterations:           ", Optim.iterations(optim_result))
    println(io, "- Run time (s):                  ", Optim.time_run(optim_result))
    println(io, "- Initial objective function value: ", PLACEHOLDER)
    println(io, "- Final objective function value:   ", PLACEHOLDER)
    println(io, " Metrics:")
    println(io, "     With objective function using P-T-X points:         Quality factor of the composition (Qcmp)")
    println(io, "     With objective function using functional relation:  ΔT(model-Henry05)")
    println(io, "- Initial metric: ", PLACEHOLDER)
    println(io, "- Final metric:   ", PLACEHOLDER)
    println(io, "\n")
    println(io, "-----------------------------------------------------")
    println(io, " Optimised variables:")
    for i in eachindex(opt_params)
        println(io, "    ", opt_params_names[i], " = ", opt_params[i])
    end
    println(io, "-----------------------------------------------------")
end
