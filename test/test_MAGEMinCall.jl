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