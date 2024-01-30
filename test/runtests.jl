using ThermoFit
using Test

CST = global_constants()
PARAMS = global_parameters()

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

    algorithm = "NelderMead";
    number_iterations_max = 1;               # set to 1 for test
    normalization = true;
    number_constraints_max = 10;             # maximum number of constraints to use for the inversion
    max_time_seconds = 5;                    # maximum time in seconds for the inversion

    JOB = job("mp", "bi", w_names, w_initial_values, w_lower_bounds, w_upper_bounds, algorithm, number_iterations_max, normalization, number_constraints_max, max_time_seconds)


    job_check_consistency(JOB)

    variables_optim, variables_optim_bounds, variables_optim_coordinates = get_variables_optim(JOB)

    w_g = calculate_w_g(variables_optim,variables_optim_coordinates, 8, 700, JOB)

    println(w_g[3])
    println(w_g[5])

    @test w_g[3] ≈ 10 + 0.1 * (700 + 273.15) + 3 * 8
    @test w_g[5] ≈ 8 + 0.2 * (700 + 273.15) + 4 * 8

end