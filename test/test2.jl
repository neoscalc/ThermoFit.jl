using ThermoFit
using CSV
using DataFrames

CST = global_constants()
PARAM = global_parameters()

function load_constraints(path, bulk_oxides)
    load_bulk(path)

    return constraints
end


function load_mineral_composition(path, mineral_elements)
    df = DataFrame(CSV.File(path, header=1))

    # rename colums with synonyms in CST.element_synonyms if element is present
    for element in keys(CST.element_synonyms)
        if element in names(df)
            rename!(df, element => CST.element_synonyms[element])
        end
    end

    # extract matrix of mineral compositions from dataframe ordered by mineral_elements
    mineral_compositions = zeros(size(df)[1], length(mineral_elements))
    for i = eachindex(mineral_elements)
        mineral_compositions[:, i] = df[!, mineral_elements[i]]
    end

    return mineral_compositions
end


function load_bulk_composition(path, bulk_oxides)
    df = DataFrame(CSV.File(path, header=1))

    # rename colums with synonyms in CST.element_synonyms if element is present
    for element in keys(CST.element_synonyms)
        if element in names(df)
            rename!(df, element => CST.element_synonyms[element])
        end
    end

    element_idx = Array{Int64}(undef, length(names(df)))
    for i = eachindex(names(df))
        element_idx[i] = CST.element_idx[names(df)[i]]
    end

    return bulk_composition
end



# test
CST = global_constants()
PARAM = global_parameters()

path = "test/data/bulk.csv"
df = DataFrame(CSV.File(path, header=1))

if PARAM.debug
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

if PARAM.debug
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

if PARAM.debug
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