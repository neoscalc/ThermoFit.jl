using ThermoFit
using Test

@testset "ThermoCheck" begin

@testset "pixelmap.jl" begin
    temperature_vec = Vector(550.:10:580)
    pressure_vec = Vector(4.:0.5:5.5)
    comp_variables_export = ["Si", "Fe", "Mg"]
    database = "mp"
    w_g = nothing
    G_0 = [0., 0., 0., 0., 0., 0., 0.]
    phase = "bi"
    # bulk from White et al. 2014b (Fig. 8)
    bulk = [64.578, 13.651, 1.586, 5.529, 8.025, 2.943, 2.000, 0.907, 0.175, 40.]
    bulk_oxides = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "MnO", "H2O"]
    sys_in = "mol"

    min_comp, p, t = pixelmap(temperature_vec, pressure_vec, bulk, bulk_oxides, database, sys_in, comp_variables_export, phase, w_g=w_g, G_0=G_0)

    @test min_comp[1][2,2] ≈ 2.7340356340442105  atol=1e-3

    w_g = repeat([0.], 21)

    min_comp, p, t = pixelmap(temperature_vec, pressure_vec, bulk, bulk_oxides, database, sys_in, comp_variables_export, phase, w_g=w_g, G_0=G_0)

    @test min_comp[1][2,2] ≈ 2.6369670586738323 atol=1e-3

    w_g = repeat([0.], 21)
    G_0 = [2., 0., 0., 0., 0., 0., 0.]
    min_comp, p, t = pixelmap(temperature_vec, pressure_vec, bulk, bulk_oxides, database, sys_in, comp_variables_export, phase, w_g=w_g, G_0=G_0)

    @test min_comp[1][2,2] ≈ 2.6265499111832478 atol=1e-3

end

@testset "qfactor.jl" begin
    #=
    create_test_points()
    =#
    constraints = ThermoCheck.create_test_points("FPWMP")
    @test length(constraints) == 12
    @test constraints[1].pressure_GPa == 0.2
    @test constraints[1].temperature_C == 550.0
    @test constraints[4].temperature_C == 700.0
    @test constraints[7].temperature_C == 650.0
    @test constraints[5].bulk == [70.999, 0.758, 12.805, 6.342, 0.075, 3.978, 0.771, 1.481, 2.7895, 0.72933, 30.0]
    @test constraints[11].bulk_oxides == ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "O", "H2O"]

    #=
    q_asm()
    =#
    @test ThermoCheck.q_asm(["mu", "cd", "bi", "pl", "q"], ["mu", "cd", "bi", "pl", "q"]) == 100.0
    @test ThermoCheck.q_asm(["cd", "bi", "pl", "q"], ["mu", "cd", "bi", "pl", "q"]) == 100.0 * (4/5)
    @test ThermoCheck.q_asm(["mu", "cd", "bi", "pl", "q"], ["afs", "cd", "bi", "pl", "q"]) == 100.0 * (4/6)

    #=
    q_vol()
    =#
    @test ThermoCheck.q_vol(["mu", "cd", "bi", "pl", "q"], ["mu", "cd", "bi", "pl", "q"], [0.2, 0.3, 0.1, 0.35, 0.05], [0.2, 0.3, 0.1, 0.35, 0.05]) == 100.0
    @test ThermoCheck.q_vol(["mu", "cd", "bi", "pl", "q"], ["cd", "mu", "bi", "pl", "q"], [0.2, 0.3, 0.1, 0.35, 0.05], [0.3, 0.2, 0.1, 0.35, 0.05]) == 100.0
    @test ThermoCheck.q_vol(["mu", "q", "bi"], ["mu", "q"], [0.4, 0.4, 0.2], [0.5, 0.5]) == 100.0 * sqrt(sum([0.9/2 * (1-0.1/0.9)^2, 0.9/2 * (1-0.1/0.9)^2]))

    #=
    probe_test_points()
    =#
    asm_vec, modes_vec = ThermoCheck.probe_test_points(constraints, "mp", "bi")
    println("Probe test points: ", asm_vec)
    @test length(intersect(asm_vec[5], ["q", "mu", "chl", "bi", "pl", "mt", "ilm", "g", "H2O"])) == length(asm_vec[5])
    @test length(setdiff(asm_vec[5], ["q", "mu", "chl", "bi", "pl", "mt", "ilm", "g", "H2O"])) == 0
    @test modes_vec[5] ≈ [0.0933726, 0.04168995, 0.08495370, 0.010032855, 0.30101063, 0.004496160, 0.0008305619, 0.32724574, 0.13636773] atol = 1e-4

    #=
    asm_diff_vector()
    =#
    asm_vec = [["mu", "cd", "bi", "pl", "q"], ["mu", "cd", "bi", "pl", "q"]]
    asm_ref = [["mu", "cd", "bi", "pl", "q"], ["mu", "cd", "bi", "pl", "q"]]
    asm_diff_vec = ThermoCheck.asm_diff_vector(asm_vec, asm_ref)
    @test asm_diff_vec == ["no diff", "no diff"]

    asm_vec = [["mu", "cd", "bi", "pl", "q"], ["mu", "cd", "bi", "pl", "q"]]
    asm_ref = [["mu", "bi", "pl", "q"], ["afs", "cd", "bi", "pl", "q"]]
    asm_diff_vec = ThermoCheck.asm_diff_vector(asm_vec, asm_ref)
    @test asm_diff_vec == ["+cd", "+mu, -afs"]

end
end