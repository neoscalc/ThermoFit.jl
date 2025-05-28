

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
    