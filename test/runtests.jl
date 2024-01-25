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

