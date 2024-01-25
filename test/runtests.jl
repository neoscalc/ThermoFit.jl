using ThermoFit
using Test

CST = global_constants()
PARAMS = global_parameters()

@testset "utils" begin
    path = "data/biotite.csv"

    element_list, mineral_element_moles = load_mineral_composition(path)
    println("test load mineral: element_list & mineral composition in moles")
    println(element_list)
    println(mineral_element_moles[2, :])

    @test mineral_element_moles[2, :] ≈ [1.200000000000000000e+01,2.614040999999999837e+00,7.367899999999999450e-02,1.628751000000000060e+00,1.513754000000000044e+00,1.154457000000000066e+00,1.531899999999999928e-02,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,1.852643000000000040e+00]

end

@testset "Bingo" begin
    
    obs_comp = [1,1.2,3.1]
    obs_unc = [0.1,0.023,0.3]
    mod_comp = [0.9,1.3,3.4]

    qcmp = bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)

    @test qcmp ≈ 69.277480978266283 atol=1e-6

end

