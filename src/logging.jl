
"""
    print_job(job)
Prints the job parameters.
"""
function print_job(job; io=stdout)
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


function print_constraints( ; io=stdout)
    println(io, "-----------------------------------------------------")
    println(io, "          CONSTRAINTS FOR THE INVERSION")
    println(io, "-----------------------------------------------------")
    println(io, " - Number of constraints: ", )
    println(io, " - SOME INFORMATION ABOUT THEM....: ")
    println(io, "-----------------------------------------------------")
end


function print_results(; io=stdout)
    date = Dates.format(Dates.now(), "dd.mm.yyyy HH:MM:SS")

    println(io, "-----------------------------------------------------")
    println(io, "             RESULTS OF THE INVERSION")
    println(io, "-----------------------------------------------------")
    println(io, " Run terminated on: ", date)
    println(io, "-----------------------------------------------------")
    println(io, " INVERSION RUN RESULTS:")
    println(io, "- Algorithm:                      ", )
    println(io, "- Number of iterations:           ", )
    println(io, "- Run time (s):                  ", )
    println(io, "- Stopped by:                   ", )
    println(io, "- Initial objective function value: ", )
    println(io, "- Final objective function value: ", )
    println(io, "- Initial metric: ", )
    println(io, "- Final metric: ", )
    println(io, "\n")
    println(io, "-----------------------------------------------------")
    println(io, " Optimised variables:")
    println(io, " PRINT THE OPTIMISED VARIABLES HERE")
    println(io, "-----------------------------------------------------")
end
