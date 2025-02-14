#=
Here a series of pissoble loss functions should be defined.

All loss functions should be of the form:

function loss_function(x::Array{Float64,1}, y::Array{Float64,1})::Float64
    return loss::Float64
end

where x is the predicted value and y is the true value.


Couple of potentially useful losses...
for ASSEMBLAGE constraints
(1) Cross entropy

for COMPOSITION constraints
(2) q_cmp
(3) chi squared

=#


"""
    bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)
Calculate the quality factor for each given component of the selected phase. Reference: Duesterhoeft & Lanari (2020) Journal of Metamorphic Geology, 38, 527-551.

## Arguments
    - mod_comp ::AbstractVector: model composition in apfu
    - obs_comp ::AbstractVector: observed composition in apfu
    - obs_unc  ::AbstractVector: uncertainty in observed composition in apfu

## Returns
    - qcmp_phase::Float64: quality factor for the given phase
"""
function bingo_calculate_qcmp_phase(mod_comp,obs_comp,obs_unc)
    # enable print statements for debugging
    PARAMS_debug::Bool = false

    # Set parameters
    FACTOR_1 = 1.00
    FACTOR_2 = 5.00

    if PARAMS_debug
        println("\n--> bingo_calculate_qcmp_phase")
        println("FACTOR_1: ", FACTOR_1)
        println("FACTOR_2: ", FACTOR_2)
        println("mod_comp: ", mod_comp)
        println("obs_comp: ", obs_comp)
        println("obs_unc(changed): ", obs_unc)
    end

    diff_abs = abs.(mod_comp .- obs_comp)
    diff_2 = diff_abs .- obs_unc/FACTOR_1

    if PARAMS_debug
        println("diff_abs: ", diff_abs)
        println("diff_2: ", diff_2)
    end

    for i = eachindex(diff_2)
        if diff_2[i] <= 0
            diff_2[i] = 0
        elseif diff_2[i] > FACTOR_2 * obs_unc[i]
            diff_2[i] = FACTOR_2 * obs_unc[i]
        end
    end

    # Calculate Q_cmp for each element
    qcmp_small = (1 .- diff_2 ./ (FACTOR_2 .* obs_unc)) .^ (mod_comp .+ 1)

    if PARAMS_debug
        println("qcmp_small: ", qcmp_small)
    end

    qcmp_phase = sum(qcmp_small)/length(qcmp_small) * 100

    if PARAMS_debug
        println("qcmp_phase: ", qcmp_phase)
    end

    return qcmp_phase
end


"""
    bingo_generate_fake_uncertainties(obs_comp)
Generate fake mineral composition uncertainties for Bingo's Qcmp (for testing purposes).
"""
function bingo_generate_fake_uncertainties(obs_comp)
    # Set parameters
    MIN_UNCERTAINTY = 0.01
    MAX_UNCERTAINTY = 0.10

    obs_unc = 0.01057612353179061 .+ 0.004459166905399905 .* obs_comp
    obs_unc[obs_unc .< MIN_UNCERTAINTY] .= MIN_UNCERTAINTY
    obs_unc[obs_unc .> MAX_UNCERTAINTY] .= MAX_UNCERTAINTY

    return obs_unc
end


"""
    chi_squared(y::AbstractVector, y_ref::AbstractVector)::Float64

## Arguments
    - y::AbstractVector: The predicted composition
    - y_ref::AbstractVector: The composition of the constraint
"""
function chi_squared(y::AbstractVector, y_ref::AbstractVector)::Float64
    chi_squared = sum((y .- y_ref).^2 ./ y_ref)
    return chi_squared
end
