"""
    global_params(solid_solution, thermodynamic_database)

solid_solution: String
    Name of the solid solution to be used in the inversion. Must be
    a valid solid solution in the thermodynamic database.

thermodynamic_database: String
    Name of the thermodynamic database to be used in the inversion.



Sets up global parameters for the inversion problem
"""

using Parameters

@with_kw struct global_params
    solid_solution::String
    thermodynamic_database::String = "mp"
end

"""
    global_constants()

element_index: Dict{String, Int64}
    Dictionary that maps element names to their index index
    in a composition vector.

element_index_reversed: Dict{Int64, String}
    Dictionary that maps element indices to their name.

oxide_index: Dict{String, Int64}
    Dictionary that maps oxide names to their index index
    in a composition vector.

oxide_index_reversed: Dict{Int64, String}
    Dictionary that maps oxide indices to their name.

Sets up global global_constants for the inversion problem
"""
@with_kw struct global_constants
    element_index::Dict{String, Int64} = Dict([("Si" => 1),
                                               ("Ti" => 2),
                                               ("Al" => 3),
                                               ("Cr" => 4),
                                               ("Fe" => 5),
                                               ("Mn" => 6),
                                               ("Mg" => 7),
                                               ("Ca" => 8),
                                               ("Na" => 9),
                                               ("K" => 10),
                                               ("O" => 11),
                                               ("H" => 12),
                                               ("C" => 13),
                                               ("S" => 14),
                                               ("Cl" => 15)])

    element_index_reversred::Dict{Int64, String} = Dict([(1 => "Si"),
                                                         (2 => "Ti"),
                                                         (3 => "Al"),
                                                         (4 => "Cr"),
                                                         (5 => "Fe"),
                                                         (6 => "Mn"),
                                                         (7 => "Mg"),
                                                         (8 => "Ca"),
                                                         (9 => "Na"),
                                                         (10 => "K"),
                                                         (11 => "O"),
                                                         (12 => "H"),
                                                         (13 => "C"),
                                                         (14 => "S"),
                                                         (15 => "Cl")])

    oxide_index::Dict{String, Int64} = Dict([("SiO2" => 1),
                                             ("TiO2" => 2),
                                             ("Al2O3" => 3),
                                             ("Cr2O3" => 4),
                                             ("FeO" => 5),
                                             ("MnO" => 6),
                                             ("MgO" => 7),
                                             ("CaO" => 8),
                                             ("Na2O" => 9),
                                             ("K2O" => 10),
                                             ("O2" => 11),
                                             ("H2O" => 12),
                                             ("CO2" => 13),
                                             ("SO2" => 14),
                                             ("Cl2" => 15)])

    oxide_index_reversed::Dict{Int64, String} = Dict([(1 => "SiO2"),
                                                      (2 => "TiO2"),
                                                      (3 => "Al2O3"),
                                                      (4 => "Cr2O3"),
                                                      (5 => "FeO"),
                                                      (6 => "MnO"),
                                                      (7 => "MgO"),
                                                      (8 => "CaO"),
                                                      (9 => "Na2O"),
                                                      (10 => "K2O"),
                                                      (11 => "O2"),
                                                      (12 => "H2O"),
                                                      (13 => "CO2"),
                                                      (14 => "SO2"),
                                                      (15 => "Cl2")])

    bingo_fac_1::Float64 = 1.00
    bingo_fac_2::Float64 = 5.00

 end