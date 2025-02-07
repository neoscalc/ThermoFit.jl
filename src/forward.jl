"""
    forward.jl

This modules contains functions handling a forward pass trough MAGEMin's G-minimiser.

## Functions
- `forward_call()`
- `calculate_w_g()`
"""


"""
    forward_call(phase, database, constraint, w_g, sys_in, gv, z_b, DB, splx_data)

Performs a forward calculation (Gibbs energy minimisation) using MAGEMin.

## Arguments
- `phase`: The name of the phase (solid solution) whichs Margules are changed.
- `database`: The database containing the thermodynamic data.
- `constraint`: The constraint object specifying the pressure, temperature, bulk composition, bulk oxides list, mineral composition and element list for the mineral composition.
- `w_g`: A vector of Margules (W-G) which is used for minimisation in the forward call.
- `sys_in`: String specifying the input format to MAGEMin ["mol"/"wt"].
- `gv`, `z_b`, `DB`, `splx_data`: The MAGEMin variables with the pointers.

## Returns
- `out`: The output containing the results of the calculation.
"""
function forward_call(phase, database, constraint, w_g, sys_in, gv, z_b, DB, splx_data)
    pressure_kbar = constraint.pressure_GPa * 10
    temperature_C = constraint.temperature_C
    bulk = constraint.bulk
    bulk_oxides = constraint.bulk_oxides

    gv = define_bulk_rock(gv, bulk, bulk_oxides, sys_in, database)
    gv, z_b, DB, splx_data = pwm_init(pressure_kbar, temperature_C, gv, z_b, DB, splx_data)

    ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss))
    ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss)

    ss_idx = findfirst(x->x==phase, ss_names)

    unsafe_copyto!(ss_struct[ss_idx].W, pointer(w_g), ss_struct[ss_idx].n_w)

    out = pwm_run(gv, z_b, DB, splx_data)
    return out
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