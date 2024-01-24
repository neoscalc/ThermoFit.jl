using ThermoFit
using CSV
using DataFrames

CST = global_constants()

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

path = "test/data/bulk.csv"
df = DataFrame(CSV.File(path, header=1))

# rename colums with synonyms in CST.element_synonyms if element is present
for element in keys(CST.element_synonyms)
    if element in names(df)
        rename!(df, element => CST.element_synonyms[element])
    end
end

element_idx = Array{Int64}(undef, length(names(df)))
for i = eachindex(names(df))
    element_idx[i] = CST.element_index[names(df)[i]]
end

oxides = Array{String}(undef, length(names(df)))
for i = eachindex(element_idx)
    oxides[i] = CST.oxide_index_reversed[element_idx[i]]
end

# rename columns with oxides
for i = eachindex(oxides)
    rename!(df, names(df)[i] => oxides[i])
end

# calculate oxides moles by dividing element moles by oxide cation number (CST.oxides_nb_cations)
for i = eachindex(names(df))
    df[!, names(df)[i]] = df[!, names(df)[i]] ./ CST.oxides_nb_cations[findfirst(x->x==names(df)[i], oxides)]

end




CST.oxides_nb_cations



mineral_elements = ["Si", "Al", "Ca", "Mg", "Fe", "K", "Na", "Ti", "O", "Mn", "H"]

load_mineral_composition("test/data/biotite.csv", mineral_elements)