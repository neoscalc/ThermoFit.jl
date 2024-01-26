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

@testset "Generation of w_g from variables_optim" begin
    w_names =  ["W(phl,annm)",
        "W(phl,obi)",
        "W(phl,east)",
        "W(phl,tbi)",
        "W(phl,fbi)",
        "W(phl,mmbi)",
        "W(annm,obi)",
        "W(annm,east)",
        "W(annm,tbi)",
        "W(annm,fbi)",
        "W(annm,mmbi)",
        "W(obi,east)",
        "W(obi,tbi)",
        "W(obi,fbi)",
        "W(obi,mmbi)",
        "W(east,tbi)",
        "W(east,fbi)",
        "W(east,mmbi)",
        "W(tbi,fbi)",
        "W(tbi,mmbi)",
        "W(fbi,mmbi)"]

    w_initial_values = [12  0  0 ;
                4  0  0 ;
                10  0.1  3 ;
                30  0  0 ;
                8  0.2  4 ;
                9  0  0 ;
                8  0  0 ;
                15  0  0 ;
                32  0  0 ;
                13.6  0  0 ;
                6.3  0  0 ;
                7  0  0 ;
                24  0  0 ;
                5.6  0  0 ;
                8.1  0  0 ;
                40  0  0 ;
                1  0  0 ;
                13  0  0 ;
                40  0  0 ;
                30  0  0 ;
                11.6  0  0]

    w_lower_bounds =   [0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0]

    w_upper_bounds =   [0 0 0;
                0 0 0;
                60 1 20;
                0 0 0;
                60 1 20;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0;
                0 0 0]

    JOB = job("mp", "bi", w_names, w_initial_values, w_lower_bounds, w_upper_bounds)

    job_check_consistency(JOB)

    variables_optim, variables_optim_bounds, variables_optim_coordinates = get_variables_optim(JOB)

    w_g = calculate_w_g(variables_optim,variables_optim_coordinates, 8, 700, JOB)

    @test w_g[3] ≈ 10 + 0.1 * 700 + 3 * 8
    @test w_g[5] ≈ 8 + 0.2 * 700 + 4 * 8

end