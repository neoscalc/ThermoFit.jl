"""
    forward.jl

This modules contains functions handling a forward pass trough MAGEMin's G-minimiser.
"""

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


