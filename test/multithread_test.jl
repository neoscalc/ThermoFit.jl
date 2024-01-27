using Base.Threads
using MAGEMin_C

MAGEMin_db = Initialize_MAGEMin("mp")


bulk = [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0]
bulk_oxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"]

sys_in = "wt"
pressures = 9 .* ones(3)
temperatures = 650 .* ones(3)

phase = "bi"
w_g_vec = [[12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6],
       [12.0, 4.0, 50.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6],
       [12.0, 4.0, 90.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]]

w_g_single = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

###
# From: https://github.com/ComputationalThermodynamics/MAGEMin/blob/main/julia/MAGEMin_wrappers.jl#L288

    # Currently, there seem to be some type instabilities or something else so that
    # some compilation happens in the threaded loop below. This interferes badly
    # in some weird way with (libsc, p4est, t8code) - in particular on Linux where
    # we get segfaults. To avoid this, we force serial compilation by calling MAGEMin
    # once before the loop.

# test runs without, but re-activate if run into SegFaults
###
# let id      = 1
#     gv          = MAGEMin_db.gv[id]
#     z_b         = MAGEMin_db.z_b[id]
#     DB          = MAGEMin_db.DB[id]
#     splx_data   = MAGEMin_db.splx_data[id]

#     point_wise_minimization(pressures[1], temperatures[1], gv, z_b, DB, splx_data, sys_in)
# end


@threads :static for i in eachindex(pressures)
    id          = Threads.threadid()

    # w_g = w_g_vec[i]

    w_g = copy(w_g_single)
    w_g[3] = w_g_single[3] + (i-1) * 40
    println("w_g = ", w_g)

    MAGEMin_db.gv[id] = define_bulk_rock(MAGEMin_db.gv[i], bulk, bulk_oxides, sys_in, MAGEMin_db.db)

    gv          = MAGEMin_db.gv[id]
    z_b         = MAGEMin_db.z_b[id]
    DB          = MAGEMin_db.DB[id]
    splx_data   = MAGEMin_db.splx_data[id]

    gv, z_b, DB, splx_data = pwm_init(pressures[i], temperatures[i], gv, z_b, DB, splx_data)

    ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
    ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    ss_idx = findfirst(x->x==phase, ss_names);

    # W = unsafe_wrap(Vector{Cdouble},ss_struct[ss_idx].W, ss_struct[ss_idx].n_w);
    unsafe_copyto!(ss_struct[ss_idx].W, pointer(w_g), ss_struct[ss_idx].n_w)
    # W .= w_g

    out = pwm_run(gv, z_b, DB, splx_data);
end

Finalize_MAGEMin(MAGEMin_db)

