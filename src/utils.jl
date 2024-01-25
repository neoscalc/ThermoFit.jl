using ThermoFit

CST = global_constants()
PARAMS = global_parameters()

function load_constraints(path_bulk, path_mineral, path_pt)
    oxide_list, bulk_oxide_mole = load_bulk_composition(path_bulk)
    element_list, mineral_element_mole = load_mineral_composition(path_mineral)
    pressure_kbar, temperature_celsius = load_pt(path_pt)

    if PARAMS.debug
        println("Bulk, mineral and pt succesfully loaded.")
    end

    # check if bulk_oxide_mole, mineral_element_mole, pressure_kbar and temperature_celsius have the same n_rows
    if size(bulk_oxide_mole)[1] != size(mineral_element_mole)[1] || size(bulk_oxide_mole)[1] != length(pressure_kbar) || size(bulk_oxide_mole)[1] != length(temperature_celsius)
        println("Error: bulk_oxide_mole, mineral_element_mole, pressure_kbar and temperature_celsius have not the same n_rows")
        return
    end

    # create an empty array constraints with size(bulk_oxide_mole)[1] rows and 1 column
    constraints = Vector{constraint}(undef, size(bulk_oxide_mole)[1])
    for i = eachindex(bulk_oxide_mole[:,1])
        constraints[i] = constraint(pressure_kbar[i], temperature_celsius[i], bulk_oxide_mole[i, :], oxide_list, mineral_element_mole[i, :], element_list)
    end

    if PARAMS.debug
        println(constraints[1:3])
    end

    return constraints
end


function load_pt(path)
    df = DataFrame(CSV.File(path, header=1))

    pressure_bar = df[:,1]
    pressure_kbar = df[:,1] ./ 10^3

    temperature_celsius = df[:,2]

    return pressure_kbar, temperature_celsius
end


function load_mineral_composition(path)
    df = DataFrame(CSV.File(path, header=1))

    if PARAMS.debug
        println(names(df))
    end

    # search column index in element_definition_uppercase and replace by element_definition
    idx_element = Array{Int64}(undef, length(names(df)))
    idx_element_dataframe = Array{Int64}(undef, length(names(df)))
    for i = eachindex(names(df))
        idx_element[i] = findfirst(x->x==names(df)[i], CST.element_definition_uppercase)
        idx_element_dataframe[i] = i
        if idx_element[i] != nothing
            rename!(df, names(df)[i] => CST.element_definition[idx_element[i]])
        end
    end

    if PARAMS.debug
        println(names(df))
        println(idx_element)
    end

    element_list = names(df)

    # extract matrix of mineral compositions from DataFrame
    mineral_element_mole = Matrix(df[:, idx_element_dataframe])

    if PARAMS.debug
        println(mineral_element_mole[1:3,:])
        println(element_list)
    end

    return element_list, mineral_element_mole
end


function load_bulk_composition(path)
    df = DataFrame(CSV.File(path, header=1))

    if PARAMS.debug
        println(names(df))
    end

    # search column index in element_definition_uppercase and replace by element_definition
    idx_element = Array{Int64}(undef, length(names(df)))
    idx_element_dataframe = Array{Int64}(undef, length(names(df)))
    for i = eachindex(names(df))
        idx_element[i] = findfirst(x->x==names(df)[i], CST.element_definition_uppercase)
        idx_element_dataframe[i] = i
        if idx_element[i] != nothing
            rename!(df, names(df)[i] => CST.element_definition[idx_element[i]])
        end
    end

    if PARAMS.debug
        println(names(df))
        println(idx_element)
    end

    # Exlude oxygen
    element_list = names(df)
    idx_oxygen = findfirst(x->x=="O", element_list)
    element_list_no_oxygen = copy(element_list)
    deleteat!(element_list_no_oxygen, idx_oxygen)
    idx_element_no_oxygen = copy(idx_element)
    deleteat!(idx_element_no_oxygen, idx_oxygen)
    idx_element_dataframe_no_oxygen = copy(idx_element_dataframe)
    deleteat!(idx_element_dataframe_no_oxygen, idx_oxygen)

    if PARAMS.debug
        println(element_list_no_oxygen)
        println(idx_element_no_oxygen)
        println(idx_element_dataframe_no_oxygen)
    end

    oxide_list = CST.oxides_definition[idx_element_no_oxygen]
    oxide_list = push!(oxide_list, "O")

    # create an array bulk_oxide_mole with idx_element_dataframe_no_oxygen columns and size(df)[1] rows
    bulk_oxide_mole = zeros(size(df)[1], length(idx_element_dataframe))

    for i = 1:size(df)[1]
        bulk_oxide_mole[i,1:end-1] = Vector(df[i, idx_element_dataframe_no_oxygen]) ./ CST.oxides_nb_cations[idx_element_no_oxygen]
        bulk_nb_oxygen = bulk_oxide_mole[i,1:end-1] .* CST.oxides_nb_oxygen[idx_element_no_oxygen]
        excess_oxygen = df[i,idx_oxygen] - sum(bulk_nb_oxygen)
        bulk_oxide_mole[i,end] = excess_oxygen
    end

    # oxide_list & bulk_oxide_mole
    return oxide_list, bulk_oxide_mole
end

function calc_structural_formula_element_from_output(out,ss_name,oxygen_norm)

    if PARAMS.debug
        println("--> calc_structural_formula_element_from_output")
    end

    # Extract the oxide names from the output
    oxides = out.oxides

    # find the indices of the oxides in CTS.oxides_definition (Philip method)
    oxide_idx = Array{Int64}(undef, length(oxides))
    for i = 1:length(out.oxides)
        oxide_idx[i] = CST.oxide_index[oxides[i]]
    end
    if PARAMS.debug
        println("Oxide index: $oxide_idx")
    end

    # Extract the solid solution index
    ss_idx = findfirst(x->x==ss_name, out.ph)
    if PARAMS.debug
        println("SS index: $ss_idx")
    end

    # Extract the solid solution composition in mol oxide
    comp_oxide = out.SS_vec[ss_idx].Comp
    if PARAMS.debug
        println("Comp oxide: $comp_oxide")
    end

    # Calculate composition in mol element
    comp_elem = comp_oxide .* CST.oxides_nb_cations[oxide_idx]
    if PARAMS.debug
        println("Comp elem: $comp_elem")
    end

    # Calculate the equivalent mol oxygen
    comp_oxygen = comp_oxide .* CST.oxides_nb_oxygen[oxide_idx]
    if PARAMS.debug
        println("Comp oxygen: $comp_oxygen")
    end

    # Calculate the normalization factor and the structural formula
    norm_factor = oxygen_norm / sum(comp_oxygen)
    comp_structural_formula = comp_elem .* norm_factor
    if PARAMS.debug
        println("Comp structural formula: $comp_structural_formula")
    end

    # find "O" and "H2O" in out.oxides and exclude them from oxide_idx, oxides and comp_structural_formula
    idx_o = findfirst(x->x=="O", oxides)
    idx_h2o = findfirst(x->x=="H2O", oxides)
    oxide_idx = deleteat!(oxide_idx, [idx_o, idx_h2o])
    oxides = deleteat!(oxides, [idx_o, idx_h2o])
    comp_structural_formula = deleteat!(comp_structural_formula, [idx_o, idx_h2o])
    if PARAMS.debug
        println("Oxide index: $oxide_idx")
        println("Oxides: $oxides")
        println("Comp structural formula: $comp_structural_formula")
    end


    return comp_structural_formula, oxides
end


function fix_order_structural_formula(comp_structural_formula_clean, oxides, constraint_element)

    if PARAMS.debug
        println("--> fix_order_structural_formula")
    end

    # Find the indices of the constraint elements in the structural formula
    constraint_element_idx = Array{Int64}(undef, length(constraint_element))
    for i = 1:length(constraint_element)
        constraint_element_idx[i] = CST.element_index[constraint_element[i]]
    end
    if PARAMS.debug
        println("Constraint element index: $constraint_element_idx")
    end

    # Obtain the corresponding indices of the oxide
    oxide_ordered = Array{String}(undef, length(constraint_element))
    oxide_idx = Array{Int64}(undef, length(constraint_element))
    for i = 1:length(constraint_element)
        oxide_ordered[i] = CST.oxide_index_reversed[constraint_element_idx[i]]
        oxide_idx[i] = findfirst(x->x==oxide_ordered[i], oxides)
    end
    if PARAMS.debug
        println("Constraint oxide (calc):  $oxide_ordered")
        println("Constraint element (def): $constraint_element")
    end

    # Adjust the order of the structural formula
    comp_structural_formula_clean_ordered = comp_structural_formula_clean[oxide_idx]
    if PARAMS.debug
        println("Comp structural formula original: $comp_structural_formula_clean")
        println("Comp structural formula ordered:  $comp_structural_formula_clean_ordered")
    end

    return comp_structural_formula_clean_ordered
end


function job_check_consistency(JOB)
    println("*** Checking job consistency ***")
    nb_wg = length(JOB.w_names)
    if isequal(size(w_initial_values), (nb_wg, 3)) == false
        error("w_initial_values must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_initial_values (size): ok")
    end
    if isequal(size(w_lower_bounds), (nb_wg, 3)) == false
        error("w_lower_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_lower_bounds (size): ok")
    end
    if isequal(size(w_upper_bounds), (nb_wg, 3)) == false
        error("w_upper_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_upper_bounds (size): ok")
    end

    println(" - Variables to be optimized [name type start min max]:")
    type_w = ["WH","WS","WV"]
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                println("    ", JOB.w_names[i], "  \t", type_w[j], "\t ", w_initial_values[i,j], "\t ", w_lower_bounds[i,j], " \t ", w_upper_bounds[i,j]) 
            end   
        end
    end
end