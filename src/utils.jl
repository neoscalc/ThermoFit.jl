CST = global_constants()
PARAMS = global_parameters()

@with_kw struct constraint
    pressure::Float64
    temperature::Float64
    bulk::Vector{Float64}
    bulk_oxides::Vector{String}
    mineral_composition::Vector{Float64}
    mineral_elements::Vector{String}
end

function load_constraints(path_bulk, path_mineral, path_pt, element_list)
    oxide_list, bulk_oxide_mole = load_bulk_composition(path_bulk)
    element_list, mineral_element_mole = load_mineral_composition(path_mineral, element_list)
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


"""
    load_pt(path)

Load pressure-temperature data from a CSV file.

# Arguments
- `path`: The path to the CSV file.

# Returns
- `pressure_kbar`: Pressure in kilobars.
- `temperature_celsius`: Temperature in degree Celsius.

"""
function load_pt(path)
    df = DataFrame(CSV.File(path, header=1))

    pressure_bar = df[:,1]
    pressure_kbar = pressure_bar ./ 10^3

    temperature_celsius = df[:,2]

    return pressure_kbar, temperature_celsius
end


"""
    load_mineral_composition(path, element_list)

Load mineral composition data from a CSV file. CSV file is in "theriak"-like format: "SI", "AL", ... and must be converted to "Si", "Al", ... after loading.
Reorder the elements according to `element_list`.

# Arguments
- `path`: The path to the CSV file containing the mineral composition data.
- `element_list`: A list of elements to extract from the CSV file and use for the inversion.

# Returns
- `element_list`: The list of elements extracted from the CSV file.
- `mineral_element_moles`: A matrix representing the mineral compositions ordered according to `element_list`.

"""
function load_mineral_composition(path, element_list)
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

    if PARAMS.debug
        println(df[2,:])
    end

    # reorder elements in df according to the order in the element_list
    mineral_element_moles_reorderd = zeros(size(df)[1], length(element_list))
    for i = 1:size(df)[1]
        for j = eachindex(element_list)
            mineral_element_moles_reorderd[i,j] = df[i, element_list[j]]
        end
    end

    if PARAMS.debug
        println(mineral_element_moles_reorderd[2,:])
        println(element_list)
    end

    # extract matrix of mineral compositions from DataFrame (reordered)
    mineral_element_moles = Matrix(mineral_element_moles_reorderd)

    return element_list, mineral_element_moles
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
        println("\n--> calc_structural_formula_element_from_output")
    end

    # Extract the oxide names from the output
    oxides = out.oxides
    if PARAMS.debug
        println("Oxides (MAGEMin): $oxides")
    end

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
        println("\n--> fix_order_structural_formula")
    end

    # Find the indices of the constraint elements in the structural formula
    constraint_element_idx = Array{Int64}(undef, length(constraint_element))
    for i = eachindex(constraint_element)
        constraint_element_idx[i] = CST.element_index[constraint_element[i]]
    end
    if PARAMS.debug
        println("Constraint element index: $constraint_element_idx")
    end

    # Obtain the corresponding indices of the oxide
    oxide_ordered = Array{String}(undef, length(constraint_element))
    oxide_idx = Array{Int64}(undef, length(constraint_element))
    for i = eachindex(constraint_element)
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
    println("\n*** CHECKING JOB CONSISTENCY ***")
    nb_wg = length(JOB.w_names)
    if isequal(size(JOB.w_initial_values), (nb_wg, 3)) == false
        error("w_initial_values must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_initial_values (size): ok")
    end
    if isequal(size(JOB.w_lower_bounds), (nb_wg, 3)) == false
        error("w_lower_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_lower_bounds (size): ok")
    end
    if isequal(size(JOB.w_upper_bounds), (nb_wg, 3)) == false
        error("w_upper_bounds must be a matrix of size (nb_wg, 3)")
    else
        println(" - w_upper_bounds (size): ok")
    end
    
    println("   -----------------------------------------------------")
    println(" - Variables to be optimized [name type start min max]:")
    type_w = ["WH","WS","WV"]
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                println("    ", JOB.w_names[i], "  \t", type_w[j], "\t ",JOB. w_initial_values[i,j], "\t ", JOB.w_lower_bounds[i,j], " \t ", JOB.w_upper_bounds[i,j])
            end
        end
    end
    println("   -----------------------------------------------------")
    println(" - Maximum number of constraints: ", JOB.number_constraints_max)
end


function get_variables_optim(JOB)

    nb_wg = length(JOB.w_names)
    type_w = ["WH","WS","WV"]

    # Calculate first the number of variables
    count = 0
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                count = count + 1
            end
        end
    end

    # Initialize the variables
    variables_optim = zeros(count)
    variables_optim_bounds = zeros(count,2)
    variables_optim_coordinates = Matrix{Int64}(undef, count, 2)

    # Extract the information
    count = 0
    for i = 1:nb_wg
        for j = 1:3
            if JOB.w_upper_bounds[i,j] > JOB.w_lower_bounds[i,j]
                count = count + 1
                variables_optim[count] = JOB.w_initial_values[i,j]
                variables_optim_bounds[count,1] = JOB.w_lower_bounds[i,j]
                variables_optim_bounds[count,2] = JOB.w_upper_bounds[i,j]
                variables_optim_coordinates[count,1] = i
                variables_optim_coordinates[count,2] = j
            end
        end
    end


    return variables_optim, variables_optim_bounds, variables_optim_coordinates
end

function calculate_w_g(variables_optim,variables_optim_coordinates, P, T, JOB)
    w_all = JOB.w_initial_values

    # replace the values using the coordinates
    for i = eachindex(variables_optim)
        w_all[variables_optim_coordinates[i,1], variables_optim_coordinates[i,2]] = variables_optim[i]
    end

    # convert deg C to K
    T_K = T + 273.15

    # calculate the g values from WG = WH + T*WS + P*WV
    w_g = w_all[:,1] .+ T_K .* w_all[:,2] .+ P .* w_all[:,3]

    return w_g
end