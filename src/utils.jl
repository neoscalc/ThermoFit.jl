CST = global_constants()
PARAMS = global_parameters()

function job_check_consistency(JOB)
    println("\n        *** CHECKING JOB CONSISTENCY ***\n")
    nb_wg = length(JOB.w_names)
    if isequal(size(JOB.w_initial_values), (nb_wg, 3)) == false
        error("w_initial_values must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_initial_values (size):        ok")
    end
    if isequal(size(JOB.w_lower_bounds), (nb_wg, 3)) == false
        error("w_lower_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_lower_bounds (size):          ok")
    end
    if isequal(size(JOB.w_upper_bounds), (nb_wg, 3)) == false
        error("w_upper_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_upper_bounds (size):          ok")
    end
    
    println("   -----------------------------------------------------")
    println(" - Variables to be optimized [name type start min max]:")
    type_w = ["WH","WS","WV"]
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                println("    ", JOB.w_names[i], "  \t", type_w[j], "\t ",JOB. w_initial_values[i,j], "\t ", JOB.w_lower_bounds[i,j], " \t ", JOB.w_upper_bounds[i,j]) 
            end
        end
    end
    println("   -----------------------------------------------------")
    println(" - Algorithm:                      ", JOB.algorithm)
    println(" - Maximum number of iterations:   ", JOB.number_iterations_max)
    println(" - Normalization:                  ", JOB.normalization)
    println(" - Maximum number of constraints:  ", JOB.number_constraints_max)
end


function get_variables_optim(JOB)

    nb_wg = length(JOB.w_names)
    type_w = ["WH","WS","WV"]

    # Calculate first the number of variables
    count = 0
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                count = count + 1
            end
        end
    end

    # Initialize the variables
    variables_optim = zeros(count)
    variables_optim_bounds = zeros(count,2)
    variables_optim_coordinates = Matrix{Int64}(undef, count, 2)

    # Extract the information
    count = 0
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                count = count + 1
                variables_optim[count] = JOB.w_initial_values[i,j]
                variables_optim_bounds[count,1] = JOB.w_lower_bounds[i,j]
                variables_optim_bounds[count,2] = JOB.w_upper_bounds[i,j]
                variables_optim_coordinates[count,1] = i
                variables_optim_coordinates[count,2] = j
            end
        end
    end


    return variables_optim, variables_optim_bounds, variables_optim_coordinates
end

function calculate_w_g(variables_optim,variables_optim_coordinates, P, T, JOB)
    # copy needed for multi_threading, so every thread has its own copy of w_all
    w_all = copy(JOB.w_initial_values)

    # replace the values using the coordinates
    for i = eachindex(variables_optim)
        w_all[variables_optim_coordinates[i,1], variables_optim_coordinates[i,2]] = variables_optim[i]
    end

    # convert deg C to K
    T_K = T + 273.15

    # calculate the g values from WG = WH + T*WS + P*WV
    w_g = w_all[:,1] .+ T_K .* w_all[:,2] .+ P .* w_all[:,3]

    return w_g
end
