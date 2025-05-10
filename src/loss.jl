#=
Here a series of possible loss functions should be defined.

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


function loss_Qfactor(mineral_composition             ::Vector{<:AbstractFloat},
                      mineral_composition_model       ::Vector{<:AbstractFloat};
                      σ_mineral_composition           ::Union{Nothing, Vector{<:AbstractFloat}} = nothing,
                      FACTOR_1                        ::Real = 1.0,
                      FACTOR_2                        ::Real = 6.0)
    
    loss_Qcmp = 100 .- quality_factor(mineral_composition,
                                      mineral_composition_model,
                                      σ_mineral_composition = σ_mineral_composition,
                                      FACTOR_1 = FACTOR_1,
                                      FACTOR_2 = FACTOR_2)
    return loss_Qcmp
end


"""
Quality factor (Bingo's Qcmp) for each given component of the selected phase.

If not measured, this implementation uses an approximation of σ(mineral composition) using "fake" uncertainties.
    
After [Duesterhoeft & Lanari (2020)](https://doi.org/10.1111/jmg.12538)
"""
function quality_factor(mineral_composition             ::Vector{<:AbstractFloat},
                        mineral_composition_model       ::Vector{<:AbstractFloat};
                        σ_mineral_composition           ::Union{Nothing, Vector{<:AbstractFloat}} = nothing,
                        FACTOR_1                        ::Real = 1.0,
                        FACTOR_2                        ::Real = 6.0)

    if isnothing(σ_mineral_composition)
        # generate fake uncertainties
        σij = bingo_generate_fake_uncertainties(mineral_composition)
    else
        σij = σ_mineral_composition
    end

    Δij = abs.(mineral_composition_model .- mineral_composition)

    Qij = zeros(length(Δij))

    for i in eachindex(Δij)
        # Let Q factor plateau within FACTOR_1 .* σij at 100.0
        if Δij[i] <= FACTOR_1 .* σij[i]
            Qij[i] = 100.
        elseif Δij[i] >= FACTOR_2 .* σij[i]
            Qij[i] = 0.
        else
            Qij[i] = 100. * (1 - (Δij[i] - FACTOR_1 * σij[i])/(FACTOR_2 * σij[i]))^(mineral_composition_model[i] + 1)
        end
    end

    Qcmp = sum(Qij)/length(Qij)

    return Qcmp
end


"""
Generate fake mineral composition uncertainties for the Quality factor (for testing purposes).
"""
function bingo_generate_fake_uncertainties(mineral_composition)
    # Constants to keep fake uncertainties within "reasonable limits"
    MIN_UNCERTAINTY = 0.01
    MAX_UNCERTAINTY = 0.10

    obs_unc = 0.01057612353179061 .+ 0.004459166905399905 .* mineral_composition
    obs_unc[obs_unc .< MIN_UNCERTAINTY] .= MIN_UNCERTAINTY
    obs_unc[obs_unc .> MAX_UNCERTAINTY] .= MAX_UNCERTAINTY

    return obs_unc
end


"""
    chi_squared(y::AbstractVector, y_ref::AbstractVector)::Float64

## Arguments
    - y       ::AbstractVector: The predicted composition
    - y_ref   ::AbstractVector: The composition of the constraint
"""
function chi_squared(y::AbstractVector, y_ref::AbstractVector)::Float64
    chi_squared_vec = (y .- y_ref).^2 ./ y_ref
    # //NOTE This is a bit a hacky fix...
    # Chi-squared loss is not defined for any entriy in y_ref = 0. This should only be a probkem when "inverting" with simulated constraints,
    # as for natural constraints compositional variables that are not measured will most likely never be considered
    # replace NaNs with 0
    chi_squared_vec[isnan.(chi_squared_vec)] .= 0

    chi_squared = sum(chi_squared_vec)
    return chi_squared
end

"""
Misfit function of the form ϕ  = (T_{model} - T_Ti-in-Bt(biotite_composition_model))^2
"""
function Ti_in_Bt_misfit(biotite_composition_model, T_C, elements)
    idx_Ti = findfirst("Ti" .== elements)
    idx_Mg = findfirst("Mg" .== elements)
    idx_Fe = findfirst("Fe" .== elements)

    Ti_apfu = biotite_composition_model[idx_Ti]
    Mg_apfu = biotite_composition_model[idx_Mg]
    Fe_apfu = biotite_composition_model[idx_Fe]

    ϕ = abs(T_C .- Ti_in_Bt_Henry05(Ti_apfu, Mg_apfu, Fe_apfu))

    return ϕ
end

"""
Ti-in-Bt thermometer of [Henry et al. (2005)](https://doi.org/10.2138/am.2005.1498)
"""
function Ti_in_Bt_Henry05(Ti_apfu, Mg_apfu, Fe_apfu)
    # function parameters of Henry at al.'s fit:
    a = -2.354
    b = 4.6482e-9
    c = -1.7283

    # convert the Ti[apfu] in 11-Ox. base to 22-Ox. base
    Ti_22Ox_base = 2 .* Ti_apfu
    # calculate XMg as (Mg)/(Mg + Fetot)
    X_Mg = Mg_apfu ./ (Mg_apfu .+ Fe_apfu)

    T_C_cubed = ((log.(Ti_22Ox_base) .- a .- c .* X_Mg.^3) ./ b)
    # if T_C_cubed is <0, Henry's thermometer is not defined, and the temperature is reported as missing
    # NOTE - This might cause a problem!
    if ismissing(T_C_cubed)
        T_C = missing
    elseif T_C_cubed < 0
        T_C = missing
    else
        T_C = T_C_cubed .^(1/3)
    end

    return T_C
end

"""
Ti saturation in biotite after [Henry et al. (2005)](https://doi.org/10.2138/am.2005.1498)

Adjusted for 11-Ox. base biotite, instead of 22-Ox. --> multiply by 0.5
"""
function Ti_saturation_Henry05(T_C, Mg_apfu, Fe_apfu)
    # function parameters of Henry at al.'s fit:
    a = -2.354
    b = 4.6482e-9
    c = -1.7283
    
    # calculate XMg as (Mg)/(Mg + Fetot)
    X_Mg = Mg_apfu ./ (Mg_apfu .+ Fe_apfu)

    return 0.5 * exp(a + b * T_C^3 + c * X_Mg^3)
end


function Ti_sat_misfit(biotite_composition_model, T_C, assemblage, elements)
    idx_Ti = findfirst("Ti" .== elements)
    idx_Mg = findfirst("Mg" .== elements)
    idx_Fe = findfirst("Fe" .== elements)

    Ti_apfu = biotite_composition_model[idx_Ti]
    Mg_apfu = biotite_composition_model[idx_Mg]
    Fe_apfu = biotite_composition_model[idx_Fe]

    if any(x -> x in assemblage, ["ilm", "ru"])
        ϕ = abs(Ti_apfu - Ti_saturation_Henry05(T_C, Mg_apfu, Fe_apfu))
    else
        Δsat = Ti_apfu - Ti_saturation_Henry05(T_C, Mg_apfu, Fe_apfu)
        if Δsat > 0
            ϕ = abs(Δsat)
        else
            ϕ = 0
        end
    end
    return ϕ   
end
