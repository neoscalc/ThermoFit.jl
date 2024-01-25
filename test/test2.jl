using ThermoFit
using CSV
using DataFrames


# test load mineral
CST = global_constants()
PARAMS = global_parameters()

path = "test/data/biotite.csv"

element_list, mineral_element_mole = load_mineral_composition(path)
println(element_list)
println(mineral_element_mole[1:5, :])

# test load bulk
CST = global_constants()
PARAMS = global_parameters()

path = "test/data/bulk.csv"

oxide_list, bulk_oxide_mole = load_bulk_composition(path)
println(oxide_list)
println(bulk_oxide_mole[1:5, :])


# test load constraints
CST = global_constants()
PARAMS = global_parameters()

path_mineral = "test/data/biotite.csv"
path_bulk = "test/data/bulk.csv"
path_pt = "test/data/pt.csv"

constraints = load_constraints(path_bulk, path_mineral, path_pt)
