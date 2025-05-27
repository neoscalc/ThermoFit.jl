
"""
    objective_function(*args)

Add an objective function Doc here.
"""
function objective_function(x0, job, constraints, nb_constraints, MAGEMin_db; loss_f::Function=loss_Qfactor, metric_f::Function=quality_factor, progress_log_io=progress_log_io)
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
    metric_vec            = zeros(nb_constraints)
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
            metric_vec[i] = 0
            loss_vec[i] = 100       #NOTE - Should be improved to consider if phase is stable and should not be (allow to pass "negative constaints"?)
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
            metric_vec[i] = metric_f(constraint_composition, composition_predicted)            
        end

    end

    # calculate the mean loss over all constraints and the fraction of constraints where the phase optimised is predicted stable
    mean_loss = sum(loss_vec) / nb_constraints
    frac_phase_present = (sum(phase_pred_stable) / nb_constraints) * 100
    # calculate the mean between the two
    residual = .5 * (mean_loss + frac_phase_present)

    # Calculate the mean Q_cmp metric
    metric = sum(metric_vec) / nb_constraints

    if !isnothing(progress_log_io)
        println(progress_log_io, "$variables_optim_local; $(residual); $(qcmp); $(frac_phase_present)")
    end

    if job.verbose
        println("\n   Residual = ", residual)
        println("   Metrics = ", metric)
        println("   Fraction of constraints where the phase is predicted stable = ", frac_phase_present)
        println("   Loss composition = ", mean_loss)
        println("\n   Optimied variables = ", variables_optim_local)
    end

    return residual
end


"""
    objective_function_func_relation(*args)

Add an objective function Doc here.
"""
function objective_function_func_relation(x0, job, constraints, nb_constraints, MAGEMin_db; loss_f::Function=Ti_sat_misfit, metric_f::Function= Ti_in_Bt_misfit, progress_log_io=progress_log_io)
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
    metric_vec          = Array{Union{Missing, Float64}}(undef, nb_constraints)

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
            loss_vec[i] = 1e20 #NOTE - Should be improved to consider if phase is stable and should not be (allow to pass "negative constaints"?)
        else
            composition_predicted = out.SS_vec[findfirst(x->x==job.phase_to_be_optimised, out.ph)].Comp_apfu

            loss_vec[i] = loss_f(composition_predicted, constraint.temperature_C, out.ph, out.elements)
            metric_vec[i] = metric_f(composition_predicted, constraint.temperature_C, out.elements)          
        end

    end

    # calculate the mean loss over all constraints and the fraction of constraints where the phase optimised is predicted stable
    residual = sum(loss_vec) / nb_constraints
    metric   = sum(skipmissing(metric_vec)) / nb_constraints

    if !isnothing(progress_log_io)
        println(progress_log_io, "$variables_optim_local; $(residual); $(metric)")
        flush(progress_log_io)
    end

    if job.verbose
        println("\n   Residual = ", residual)
        println("   ΔT(model-Henry05) = ", metric)
        println("\n   Optimied variables = ", variables_optim_local)
    end

    return residual
end
