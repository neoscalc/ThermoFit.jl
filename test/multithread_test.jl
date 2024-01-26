using Base.Threads
using MAGEMin_C

MAGEMin_db = Initialize_MAGEMin("mp")

bulk = [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0]
bulk_oxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"]

sys_in = "wt"
pressures = [8.0, 8.0]
temperatures = [650, 600]

let id      = 1
    gv          = MAGEMin_db.gv[id]
    z_b         = MAGEMin_db.z_b[id]
    DB          = MAGEMin_db.DB[id]
    splx_data   = MAGEMin_db.splx_data[id]

    point_wise_minimization(pressures[1], temperatures[1], gv, z_b, DB, splx_data, sys_in)
end


@threads :static for i in eachindex(pressures)
    id          = Threads.threadid()
    MAGEMin_db.gv[id] = define_bulk_rock(MAGEMin_db.gv[i], bulk, bulk_oxides, sys_in, MAGEMin_db.db)

    gv          = MAGEMin_db.gv[id]
    z_b         = MAGEMin_db.z_b[id]
    DB          = MAGEMin_db.DB[id]
    splx_data   = MAGEMin_db.splx_data[id]


    out = point_wise_minimization(pressures[i], temperatures[i], gv, z_b, DB, splx_data)

    # gv, z_b, DB, splx_data = pwm_init(pressure, temperature, gv, z_b, DB, splx_data);
    # out = pwm_run(gv, z_b, DB, splx_data);
end

# finalize_MAGEMin(gv,DB, z_b)

