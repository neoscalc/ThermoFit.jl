"""
    forward.jl

This modules contains functions handling a forward pass trough MAGEMin's G-minimiser.

## Functions
- `forward_call()`
- `calculate_w_g()`
"""


function prepare_forward_call(variables_optim_local, MAGEMin_db, job, constraint)
    # identify thread and acess the MAGEMin_db of the thread
    id          = Threads.threadid()
    # println("   Thread ",id," is working on constraint ",i," of ",nb_constraints)
    gv          = MAGEMin_db.gv[id]
    z_b         = MAGEMin_db.z_b[id]
    DB          = MAGEMin_db.DB[id]
    splx_data   = MAGEMin_db.splx_data[id]

    # Calculate Margules (W-G) 
    idx_margules    = findall(job.var_optim_type .== "W")
    if !isempty(idx_margules)
        w_optimised     = variables_optim_local[idx_margules]
        P_kbar          = constraint.pressure_GPa * 10
        T_K             = constraint.temperature_C + 273.15
        w_matrix        = copy(job.w_initial_values)   # copy needed for multi_threading, so every thread has its own copy of w_all
        w_coordinates   = job.margules_optim_coord
            
        w_g = calculate_w_g(w_optimised, P_kbar, T_K, w_matrix, w_coordinates)
    else
        w_g = nothing
    end

    # Calculate G0 corrections
    idx_g0_coor    = findall(job.var_optim_type .== "G0_corr")
    if !isempty(idx_g0_coor)
        g0_corr = variables_optim_local[idx_g0_coor]
        g0_corr_endmembers = job.g0_corr_endmembers
    else
        g0_corr = nothing
        g0_corr_endmembers = nothing
    end
    
    return constraint, gv, z_b, DB, splx_data, w_g, g0_corr, g0_corr_endmembers
end

"""
    forward_call(phase, database, constraint, sys_in, gv, z_b, DB, splx_data; w_g = nothing, g0_corr_vec = nothing, g0_corr_em = nothing)

Performs a forward calculation (Gibbs energy minimisation) using MAGEMin.

## Arguments
- `phase`: The name of the phase (solid solution) whichs Margules are changed.
- `database`: The database containing the thermodynamic data.
- `constraint`: The constraint object specifying the pressure, temperature, bulk composition, bulk oxides list, mineral composition and element list for the mineral composition.
- `sys_in`: String specifying the input format to MAGEMin ["mol"/"wt"].
- `gv`, `z_b`, `DB`, `splx_data`: The MAGEMin variables with the pointers.
- `w_g`: A vector of Margules (W-G) which is used for minimisation in the forward call.
- `g0_corr_vec`: A vector of corrections to the G0 values of the endmembers.
- `g0_corr_em`: A vector of endmembers corresponding to the corrections.

## Returns
- `out`: The output containing the results of the calculation.
"""
function forward_call(phase, database, constraint, gv, z_b, DB, splx_data; w_g = nothing, g0_corr_vec = nothing, g0_corr_em = nothing)
    pressure_kbar = constraint.pressure_GPa * 10
    temperature_C = constraint.temperature_C
    bulk = constraint.bulk
    bulk_oxides = constraint.bulk_oxides
    sys_in = constraint.sys_in

    gv = define_bulk_rock(gv, bulk, bulk_oxides, sys_in, database)
    gv, z_b, DB, splx_data = pwm_init(pressure_kbar, temperature_C, gv, z_b, DB, splx_data)

    ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss))
    ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss)

    ss_idx = findfirst(x->x==phase, ss_names)

    if !isnothing(w_g)
        unsafe_copyto!(ss_struct[ss_idx].W, pointer(w_g), ss_struct[ss_idx].n_w)
    end

    if !isnothing(g0_corr_vec) && !isnothing(g0_corr_em)
        ss_gbase_mod = calculate_g0_corr(g0_corr_vec, g0_corr_em, ss_struct, ss_idx)
        unsafe_copyto!(ss_struct[ss_idx].gbase, pointer(ss_gbase_mod), ss_struct[ss_idx].n_em)
    end

    out = pwm_run(gv, z_b, DB, splx_data)
    return out
end


function forward_call_multi_phase(phases            ::AbstractArray,
                                  database          ::AbstractString,
                                  constraint        ::Constraint,
                                  gv, z_b, DB, splx_data;
                                  phases_w_g        ::Union{AbstractArray,Nothing} = nothing,
                                  phases_g0_corr_vec::Union{AbstractArray,Nothing} = nothing,
                                  w_g               ::Union{AbstractArray,Nothing} = nothing,
                                  g0_corr_vec       ::Union{AbstractArray,Nothing} = nothing,
                                  g0_corr_em        ::Union{AbstractArray,Nothing} = nothing)

    # check in the pahse list if for the specific mineral a W or G0 should be altered in the call
    # then essentially move into the forward_call function

end

"""
    calculate_g0_corr(g0_corr_vec, g0_corr_em, ss_struct, ss_idx)

Adds a "correction" to the G0 values of endmembers in a solid solution.

## Arguments
- `g0_corr_vec`: A vector of 'n' corrections to be added on to G0 values of the endmembers.
- `g0_corr_em`: A corresponding vector of 'n' endmembers to which the corrections are applied.
- `ss_struct`: The solid solution structure from MAGEMin containing the endmembers.
- `ss_idx`: The index of the solid solution in the structure.
"""
function calculate_g0_corr(g0_corr_vec, g0_corr_em, ss_struct, ss_idx)
    ss_gbase = unsafe_wrap(Vector{Float64}, ss_struct[ss_idx].gbase, ss_struct[ss_idx].n_em)
    ss_emlist = unsafe_wrap(Vector{Ptr{Int8}}, ss_struct[ss_idx].EM_list, ss_struct[ss_idx].n_em)
    ss_emlist = unsafe_string.(ss_emlist)
    ss_gbase_mod = copy(ss_gbase)
    for (em, g0_corr) in zip(g0_corr_em, g0_corr_vec)
        em_idx = findfirst(ss_emlist .== em)
        ss_gbase_mod[em_idx] += g0_corr
    end
    return ss_gbase_mod
end


"""
    calculate_w_g(margules_optim, P_kbar, T_K, w_matrix, w_coordinates)

Calculates the Margules (W-G) at fixed pressure and temperature for the forward pass.
The margules variables of the inversion are inserted into the Margules matrix and W_G is calculated.

## Arguments
- `margules_optim`: The Margules variables used for the inversion.
- `P_kbar`: The pressure in kbar.
- `T_K`: The temperature in Kelvin.
- `w_matrix`: The Margules matrix (n x 3). From JOB.w_initial_values.
- `w_coordinates`: The coordinates of the Margules variables in the matrix.
"""
function calculate_w_g(margules_optim, P_kbar, T_K, w_matrix, w_coordinates)
    # replace the values using the coordinates
    w_matrix[w_coordinates] .= margules_optim

    # calculate the g values from WG = WH - T*WS + P*WV, using P[kbar] and T[K] to match MAGEMin's calculation
    w_g = w_matrix[:,1] .- T_K .* w_matrix[:,2] .+ P_kbar .* w_matrix[:,3]

    return w_g
end