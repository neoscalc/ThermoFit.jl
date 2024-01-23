"""
    global_params(solid_solution, thermodynamic_database)

solid_solution: String
    Name of the solid solution to be used in the inversion. Must be
    a valid solid solution in the thermodynamic database.

thermodynamic_database: String
    Name of the thermodynamic database to be used in the inversion.


    
Sets up global parameters for the inversion problem
"""
@with_kw struct global_params
    solid_solution::String
    thermodynamic_database::String

    

end

@with_kw struct global_constants
    element_index::Dict{String, Int64} = Dict([("Si" => 1), ("Ti" => 2)])
    # oxide_index::Dict{String, Int64}

    bingo_fac_1::Float64 = 1.00
    bingo_fac_2::Float64 = 5.00

 end