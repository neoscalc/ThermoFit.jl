"""
    forward.jl

This modules contains functions handling a forward pass trough MAGEMin's G-minimiser.
"""


"""
    forward_call(phase, database, constraint, w_g, sys_in, gv, z_b, DB, splx_data)

This function performs a forward calculation (Gibbs energy minimisation) using MAGEMin.

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
    pressure = constraint.pressure
    temperature = constraint.temperature
    bulk = constraint.bulk
    bulk_oxides = constraint.bulk_oxides

    gv = define_bulk_rock(gv, bulk, bulk_oxides, sys_in, database);
    gv, z_b, DB, splx_data = pwm_init(pressure, temperature, gv, z_b, DB, splx_data);

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

    return out
end


function callMAGEMin()

    database = "mp";                    # select database here, ig, igd, alk, mp, mb, um

    global gv, z_b, DB, splx_data   = init_MAGEMin(database);

    gv = use_predefined_bulk_rock(gv, 0, database);

    # # Print information on the selected database
    gv, z_b, DB, splx_data = pwm_init(5, 650, gv, z_b, DB, splx_data);

    ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
    ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    print("\n---------------------------------\n");
    print("    Database information (",database,")\n");
    print("----------------------------------\n");
    for i=1:gv.len_ss
        print("   ",ss_names[i],": ",i,"; n_W's, ",ss_struct[i].n_w,"; n_em's, ",ss_struct[i].n_em,"\n")
    end
    print("----------------------------------\n\n");

    # Test 1: run MAGEMin with default parameters
    gv, z_b, DB, splx_data = pwm_init(5, 650, gv, z_b, DB, splx_data);
    out       = pwm_run(gv, z_b, DB, splx_data);

    return out
end


