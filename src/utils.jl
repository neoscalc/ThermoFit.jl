using ThermoFit

CST = global_constants()

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

