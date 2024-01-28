
# function to calculate a pixel map using MAGEMin
function pixelmap_calc(ti, pi, bulk, bulk_oxides, phase, database, w_g, sys_in, variables_export)

    pressure = zeros(length(pi)*length(ti))
    temperature = zeros(length(pi)*length(ti))
    idx_p = zeros(length(pi)*length(ti))
    idx_t = zeros(length(pi)*length(ti))
    for i = 1:length(pi)
        for j = 1:length(ti)
            pressure[(i-1)*length(ti)+j] = pi[i]
            temperature[(i-1)*length(ti)+j] = ti[j]
            idx_p[(i-1)*length(ti)+j] = i
            idx_t[(i-1)*length(ti)+j] = j
        end
    end

    # Create a structure result to store the results of the minimisation
    result = DataFrame(
        pressure = pressure,
        temperature = temperature,
        idx_p = idx_p,
        idx_t = idx_t,
    )
    # Add columns in results using variables_export
    for i = 1:length(variables_export)
        result[!, variables_export[i]] = zeros(length(pi)*length(ti))
    end

    MAGEMin_db  = Initialize_MAGEMin(database, verbose=false);

    @threads for i in ProgressBar(1:length(pressure))

        # identify thread and acess the MAGEMin_db of the thread
        id          = Threads.threadid()
        # println("   Thread ",id," is working on constraint ",i," of ",nb_to_use)
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        gv = define_bulk_rock(gv, bulk, bulk_oxides, sys_in, database);
        gv, z_b, DB, splx_data = pwm_init(pressure[i], temperature[i], gv, z_b, DB, splx_data);

        ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
        ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

        ss_idx = findfirst(x->x==phase, ss_names);

        unsafe_copyto!(ss_struct[ss_idx].W, pointer(w_g), ss_struct[ss_idx].n_w)

        # println("W = ", w_g)

        # println("** calling MAGEMin at: P = $(pressure) kbar and T = $(temperature) K")
        # println("   pressure = ", pressure, " kbar")
        # println("   temperature = ", temperature, " K")
        # println("   bulk composition = ", bulk)
        # println("   bulk oxides = ", bulk_oxides)

        out = pwm_run(gv, z_b, DB, splx_data);

        if (phase in out.ph)
            comp_structural_formula_clean, oxides = calc_structural_formula_element_from_output(out,phase,12)

            # fix order of the structural formula to match variables_export
            comp_structural_formula_clean_ordered = fix_order_structural_formula(comp_structural_formula_clean, oxides, variables_export)

            # Write the results in the result structure
            for j = 1:length(variables_export)
                result[!, variables_export[j]][i] = comp_structural_formula_clean_ordered[j]
            end
        end

    end


    Finalize_MAGEMin(MAGEMin_db)

    return result

end

