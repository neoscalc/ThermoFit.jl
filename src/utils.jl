using ThermoFit

CST = global_constants()
PARAMS = global_params()

function load_bulk(path)
    df = DataFrame(CSV.File(path, header=1))

    # Create an array of the bulk composition
    bulk_comp = zeros(Float64, nrow(df), length(CST.element_index))

    O = df[:, "O"]
    Si = df[:, "SI"]
    Ti = df[:, "TI"]
    Al = df[:, "AL"] ./ 2
    Fe = df[:, "FE"]
    Mg = df[:, "MG"]
    Mn = df[:, "MN"]
    Ca = df[:, "CA"]
    Na = df[:, "NA"] ./ 2
    K = df[:, "K"] ./ 2
    H = df[:, "H"] ./ 2

    bulk_comp[:, CST.element_index["O"]] = O
    bulk_comp[:, CST.element_index["Si"]] = Si
    bulk_comp[:, CST.element_index["Ti"]] = Ti
    bulk_comp[:, CST.element_index["Al"]] = Al
    bulk_comp[:, CST.element_index["Fe"]] = Fe
    bulk_comp[:, CST.element_index["Mg"]] = Mg
    bulk_comp[:, CST.element_index["Mn"]] = Mn
    bulk_comp[:, CST.element_index["Ca"]] = Ca
    bulk_comp[:, CST.element_index["Na"]] = Na
    bulk_comp[:, CST.element_index["K"]] = K
    bulk_comp[:, CST.element_index["H"]] = H

    return bulk_comp
end

function load_mineral(path)
    df = DataFrame(CSV.File(path, header=1))

    # Create an array of the bulk composition
    bulk_comp = zeros(Float64, nrow(df), length(CST.element_index))

    O = df[:, "O"]
    Si = df[:, "SI"]
    Ti = df[:, "TI"]
    Al = df[:, "AL"]
    Fe = df[:, "FE"]
    Mg = df[:, "MG"]
    Mn = df[:, "MN"]
    Ca = df[:, "CA"]
    Na = df[:, "NA"]
    K = df[:, "K"]
    H = df[:, "H"]

    bulk_comp[:, CST.element_index["O"]] = O
    bulk_comp[:, CST.element_index["Si"]] = Si
    bulk_comp[:, CST.element_index["Ti"]] = Ti
    bulk_comp[:, CST.element_index["Al"]] = Al
    bulk_comp[:, CST.element_index["Fe"]] = Fe
    bulk_comp[:, CST.element_index["Mg"]] = Mg
    bulk_comp[:, CST.element_index["Mn"]] = Mn
    bulk_comp[:, CST.element_index["Ca"]] = Ca
    bulk_comp[:, CST.element_index["Na"]] = Na
    bulk_comp[:, CST.element_index["K"]] = K
    bulk_comp[:, CST.element_index["H"]] = H

    return biotite_comp
end

function calc_structural_formula_element_from_output(out,ss_name,oxygen_norm)
    # find the indices of the oxides in CTS.oxides_definition (Philip method)
    oxide_idx = Array{Int64}(undef, length(out.oxides))
    for i = 1:length(out.oxides)
        oxide_idx[i] = CST.oxide_index[out.oxides[i]]
    end
    if PARAMS.debug
        println("Oxide index: $oxide_idx")
    end

    ss_idx = findfirst(x->x==ss_name, out.ph)
    if PARAMS.debug
        println("SS index: $ss_idx")
    end

    comp_oxide = out.SS_vec[ss_idx].Comp
    if PARAMS.debug
        println("Comp oxide: $comp_oxide")
    end

    comp_elem = comp_oxide .* CST.oxides_nb_cations[oxide_idx]
    if PARAMS.debug
        println("Comp elem: $comp_elem")
    end
    
    comp_oxygen = comp_oxide .* CST.oxides_nb_oxygen[oxide_idx]
    if PARAMS.debug
        println("Comp oxygen: $comp_oxygen")
    end

    norm_factor = oxygen_norm / sum(comp_oxygen)
    comp_structural_formula = comp_elem .* norm_factor
    if PARAMS.debug
        println("Comp structural formula: $comp_structural_formula")
    end

    return comp_structural_formula
end