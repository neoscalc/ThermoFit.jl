

function plot_convergence(filepath          ::String,
                          convergence_log   ::String,
                          loss_f            ::String,
                          metric_f          ::String;
                          START_AT_EPOCH = 5,
                          AVERAGE_N = 5)

    # Function to compute the moving average
    moving_average(vs,n) = [mean(@view vs[i:(i+n-1)]) for i in 1:(length(vs)-(n-1))]

    if endswith(convergence_log, ".log")
        data = CSV.read(convergence_log, DataFrame, header=false, delim=';')
    elseif isdir(convergence_log)
        log_files = filter(f -> endswith(f, ".log"), readdir(convergence_log))
        data = [CSV.read(joinpath(convergence_log, f), DataFrame, header=false, delim=';') for f in log_files]
    else
        @error "Invalid convergence log file or directory: $convergence_log"
    end

    # When multiple logs (random search) find the curve with minimal loss to highlight it
    if data isa Vector{DataFrame}
        min_loss_idx = argmin([minimum(d[:, 2]) for d in data])
    end

    fig = Figure(size=(500, 600))
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])

    if data isa DataFrame
        lines!(ax1, moving_average(data[START_AT_EPOCH:end, 2], AVERAGE_N), color=:blue)
        lines!(ax2, moving_average(data[START_AT_EPOCH:end, 3], AVERAGE_N), color=:red)
    elseif data isa Vector{DataFrame}
        for df in data
            lines!(ax1, moving_average(df[START_AT_EPOCH:end, 2], AVERAGE_N), color=:blue, alpha=0.2)
            lines!(ax2, moving_average(df[START_AT_EPOCH:end, 3], AVERAGE_N), color=:red, alpha=0.2)
        end

        # Highlight the curve with minimal loss
        lines!(ax1, moving_average(data[min_loss_idx][START_AT_EPOCH:end, 2], AVERAGE_N), color=:blue, linewidth=2)
        lines!(ax2, moving_average(data[min_loss_idx][START_AT_EPOCH:end, 3], AVERAGE_N), color=:red, linewidth=2)
    end

    ax1.title = "Loss $AVERAGE_N-epoch moving average of $loss_f"
    ax1.xlabel = "Epochs"
    ax1.ylabel = "Loss value"

    ax2.title = "Metric $AVERAGE_N-epoch moving average of $metric_f"
    ax2.xlabel = "Epochs"
    ax2.ylabel = "Metric value"

    save(filepath, fig)
end


function x_optim_v_resiudal(filepath          ::String,
                            convergence_log   ::String;
                            var_names         ::Union{Nothing, Vector{<:AbstractString}} = nothing)

    if endswith(convergence_log, ".log")
        data = CSV.read(convergence_log, DataFrame, header=false, delim=';')
    elseif isdir(convergence_log)
        log_files = filter(f -> endswith(f, ".log"), readdir(convergence_log))
        data = [CSV.read(joinpath(convergence_log, f), DataFrame, header=false, delim=';') for f in log_files]
    end

    n_var_optims = length(eval(Meta.parse(data[1][1, 1])))

    # extract the optimised variables, the starting values, and the minimum residuals
    x_optim_arr = Array{Float64}(undef, 0, n_var_optims)
    x_starting_arr = Array{Float64}(undef, 0, n_var_optims)
    min_resiudal_arr = Array{Float64}(undef, length(data))
    for (i, df) in enumerate(data)
        idx_min_loss = argmin(df[:, 2])

        x_optim = eval(Meta.parse(df[idx_min_loss, 1]))
        x_optim_arr = vcat(x_optim_arr, x_optim')

        x_starting = eval(Meta.parse(df[1, 1]))
        x_starting_arr = vcat(x_starting_arr, x_starting')

        min_resiudal = df[idx_min_loss, 2]
        min_resiudal_arr[i] = min_resiudal
    end

    y_lower_lim = minimum(filter(x -> x<1e3,min_resiudal_arr))
    y_lower_lim = y_lower_lim - 0.05 * y_lower_lim
    y_upper_lim = maximum(filter(x -> x<1e3,min_resiudal_arr))
    y_upper_lim = y_upper_lim + 0.1 * y_upper_lim

    colors = cgrad(:phase, n_var_optims, categorical=true)
    x_labels = isnothing(var_names) ? ["Variable $i" for i in 1:n_var_optims] : var_names

    fig = Figure(size=(600, n_var_optims*250))

    for i in 1:n_var_optims
        ax1 = Axis(fig[i, 1])
        scatter!(ax1, x_optim_arr[:, i], min_resiudal_arr, color=colors[i], markersize=8)
        ylims!(ax1, y_lower_lim, y_upper_lim)
        ax1.xlabel = x_labels[i]
        ax1.ylabel = "Residual"

        ax2 = Axis(fig[i, 2])
        arrow_vec = hcat(x_starting_arr[:, 1], min_resiudal_arr) - hcat(x_optim_arr[:, 1], min_resiudal_arr)
        arrows!(ax2, x_starting_arr[:, i], min_resiudal_arr, arrow_vec[:, 1], arrow_vec[:, 2], color=colors[i])
        ylims!(ax2, y_lower_lim, y_upper_lim)
        ax2.xlabel = "$(x_labels[i]) starting > optimised"
        ax2.ylabel = "Residual"
    end

    save(filepath, fig)
end