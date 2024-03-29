"""
    global_params(solid_solution, thermodynamic_database)

solid_solution: String
    Name of the solid solution to be used in the inversion. Must be
    a valid solid solution in the thermodynamic database.

thermodynamic_database: String
    Name of the thermodynamic database to be used in the inversion.



Sets up global parameters for the inversion problem
"""
@with_kw struct global_parameters

    debug::Bool = false                             # Debug mode: true or false

    thermodynamic_database::String = "mp"           # Thermodynamic database to use (mp, alk, ig, igd, mb, um)
    solid_solution::String = "bi"                   # Solid solution to optimise (MAGEMin name)

end

@with_kw struct job

    thermodynamic_database::String = "mp"           # Thermodynamic database to use (mp, alk, ig, igd, mb, um)
    solid_solution::String = "bi"                   # Solid solution to optimise (MAGEMin name)

    w_names::Vector{String}                         # Names of the Margules for the solid solution
    w_initial_values::Matrix{Float64}               # Initial values of the Margules for the solid solution (n_w x 3; columns: WH, WS, WV)
    w_lower_bounds::Matrix{Float64}                 # Lower bounds of the Margules for the solid solution (n_w x 3; columns: WH, WS, WV)
    w_upper_bounds::Matrix{Float64}                 # Upper bounds of the Margules for the solid solution (n_w x 3; columns: WH, WS, WV)

    algorithm::String = "Nelder-Mead"               # Algorithm to use for the inversion (Nelder-Mead, Particle Swarm)
    number_iterations_max::Int64 = 1000             # Maximum number of iterations
    normalization::Bool = true                      # Normalize the variables to optimise (Margules) for G-minimisation
    number_constraints_max::Int64 = 100             # Number of constraints to be used during inversion (0=all)
    max_time_seconds::Int64 = 1000                  # Maximum time in seconds for the inversion (0=unlimited)

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

    element_synonyms::Dict{String, String} = Dict([("SI" => "Si"),
                                                   ("TI" => "Ti"),
                                                   ("AL" => "Al"),
                                                   ("CR" => "Cr"),
                                                   ("FE" => "Fe"),
                                                   ("MN" => "Mn"),
                                                   ("MG" => "Mg"),
                                                   ("CA" => "Ca"),
                                                   ("NA" => "Na"),
                                                   ("K" => "K"),
                                                   ("H" => "H"),
                                                   ("O" => "O"),
                                                   ("C" => "C"),
                                                   ("S" => "S"),
                                                   ("CL" => "Cl")])

    element_to_oxides::Dict{String, String} = Dict([("Si" => "SiO2"),
                                                    ("Ti" => "TiO2"),
                                                    ("Al" => "Al2O3"),
                                                    ("Cr" => "Cr2O3"),
                                                    ("Fe" => "FeO"),
                                                    ("Mn" => "MnO"),
                                                    ("Mg" => "MgO"),
                                                    ("Ca" => "CaO"),
                                                    ("Na" => "Na2O"),
                                                    ("K" => "K2O"),
                                                    ("H" => "H2O"),
                                                    ("O" => "O"),
                                                    ("C" => "CO2"),
                                                    ("S" => "SO2"),
                                                    ("Cl" => "Cl2")])

    oxide_to_elements::Dict{String, String} = Dict([("SiO2" => "Si"),
                                                            ("TiO2" => "Ti"),
                                                            ("Al2O3" => "Al"),
                                                            ("Cr2O3" => "Cr"),
                                                            ("FeO" => "Fe"),
                                                            ("MnO" => "Mn"),
                                                            ("MgO" => "Mg"),
                                                            ("CaO" => "Ca"),
                                                            ("Na2O" => "Na"),
                                                            ("K2O" => "K"),
                                                            ("H2O" => "H"),
                                                            ("O" => "O"),
                                                            ("CO2" => "C"),
                                                            ("SO2" => "S"),
                                                            ("Cl2" => "Cl")])

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
                                             ("H2O" => 11),
                                             ("O" => 12),
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
                                                      (11 => "H2O"),
                                                      (12 => "O"),
                                                      (13 => "CO2"),
                                                      (14 => "SO2"),
                                                      (15 => "Cl2")])

    element_definition::Vector{String} = ["Si", "Ti", "Al", "Cr", "Fe", "Mn", "Mg", "Ca", "Na", "K", "O", "H", "C", "S", "Cl"]
    element_definition_uppercase::Vector{String} = ["SI", "TI", "AL", "CR", "FE", "MN", "MG", "CA", "NA", "K", "O", "H", "C", "S", "CL"]
    oxides_definition::Vector{String} = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "O", "H2O", "CO2", "SO2", "Cl2"]
    oxides_nb_cations::Vector{Int64} = [1,1,2,2,1,1,1,1,2,2,0,2,1,1,2]
    oxides_nb_oxygen::Vector{Int64} = [2,2,3,3,1,1,1,1,1,1,1,1,2,2,0]


    bingo_fac_1::Float64 = 1.00
    bingo_fac_2::Float64 = 5.00

 end