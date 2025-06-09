using ThermoFit
using Test

@testset "ThermoCheck" begin

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