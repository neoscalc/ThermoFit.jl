using ThermoFit

CST = global_constants()
PARAMS = global_params()

function bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)
    # Calculate the quality factor for each given component of the selected phase
    # mod_comp: model composition in apfu
    # obs_comp: observed composition in apfu
    # obs_unc: uncertainty in observed composition in apfu
    #
    # Returns: 
    # qcmp: quality factor for the given phase
    #
    # Reference: Duesterhoeft & Lanari (2020) Journal of Metamorphic Geology, 38, 527-551.
    
    # Set parameters
    fac_1 = CST.bingo_fac_1
    fac_2 = CST.bingo_fac_2

    # Replace obs_unc values smaller than 0.01 with 0.01
    obs_unc[obs_unc .< 0.01] .= 0.01

    diff_abs = abs.(mod_comp .- obs_comp)
    diff_2 = diff_abs .- obs_unc/fac_1

    for i = 1:length(diff_2)
        if diff_2[i] <= 0
            diff_2[i] = 0
        elseif diff_2[i] > fac_2 * obs_unc[i]
            diff_2[i] = fac_2 * obs_unc[i]
        end
    end

    # Calculate Q_cpm for each element
    qcmp_small = (1 .- diff_2 ./ (fac_2 .* obs_unc)) .^ (mod_comp .+ 1)

    print("qcmp_small:")
    println(qcmp_small)

    qcmp_phase = sum(qcmp_small)/length(qcmp_small) * 100
    
    return qcmp_phase
end

function bingo_generate_fake_uncertainties(obs_comp)
    # Generate fake mineral composition uncertainties for Bingo
    #
    # Returns:
    # obs_unc: fake uncertainties for testing
    #
    # Reference: 

    # Set parameters
    min_uncertainty = 0.01
    max_uncertainty = 0.10

    obs_unc = 0.01057612353179061 .+ 0.004459166905399905 .* obs_comp
    obs_unc[obs_unc .< min_uncertainty] .= min_uncertainty
    obs_unc[obs_unc .> max_uncertainty] .= max_uncertainty

    return obs_unc
end


# The code below is used to generate the parameters for the fake uncertainties
test_comp = [0.01,0.37,3.27,6.50,20.6]
test_unc =  [0.01,0.01,0.02,0.05,0.10]

test_unc./test_comp

using Plots
plot(test_comp,test_unc,seriestype=:scatter)

# fit data (to be replaced b sigmoid function?)
using Polynomials
p = fit(test_comp, test_unc, 1)
plot!([0,25], p.([0,25]), label = "fit", color=:red, linewidth=2)

p
