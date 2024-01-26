using ThermoFit
callMAGEMin()

# test objective_function
database = "mp"

constraint_A = Constraint(8.0, 650, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])
constraint_B = Constraint(8.0, 600, [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0], ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"])

constraints = [constraint_A, constraint_B]

phase = "bi"
margules = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

@btime objective_function(margules, constraints, database, phase, "wt")
# 213.704 ms (896 allocations: 58.55 KiB)



# ----------------------------------------------------------------------------------------------------------------------------
using ThermoFit

global gv, z_b, DB, splx_data  = init_MAGEMin("mp");

pressure = 6.997
temperature = 533.0
bulk = [17.70151, 0.017186, 1.903342, 0.116005, 0.079485, 0.002148, 0.272827, 0.693882, 0.691734, 10.741208, 0.00580100000001238]
# bulk = [17.507749, 0.167715, 1.7533885, 1.049855, 0.361568, 0.008712, 0.108906, 0.2940465, 0.514037, 10.8906125, 0.052494500000001665]
bulk_oxides = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "MnO", "CaO", "Na2O", "K2O", "H2O", "O"]

gv = define_bulk_rock(gv, bulk, bulk_oxides, "wt", "mp");

gv.verbose  = 1
# gv.solver = 0           # Legacy, 0. PGE, 1 (default)

gv, z_b, DB, splx_data = pwm_init(pressure, temperature, gv, z_b, DB, splx_data);

out = pwm_run(gv, z_b, DB, splx_data);

finalize_MAGEMin(gv,DB, z_b)



# Try pass gamma: [-944.220506,-1740.794181,-777.187108,-672.871086,-343.625254,-876.509241,-806.168923,-989.746167,-226.716701,-491.516729,-347.998070]
