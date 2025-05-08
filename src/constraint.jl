struct Constraint
    pressure_GPa                    ::Float64
    temperature_C                   ::Float64
    bulk                            ::Vector{Float64}
    bulk_oxides                     ::Vector{String}
    sys_in                          ::String
    assemblage                      ::Union{Nothing, Vector{String}}
    mineral_composition_apfu        ::Union{Nothing, OrderedDict{String, Vector{Float64}}}
    mineral_elements                ::Vector{String}
end


"""
    read_constraints_from_yml(path::AbstractString)

Reads constraints from a YAML file.

## Arguments
- `path`: The path to the YAML file.

## Returns
- `constraints_vec`: A vector of `Constraint` objects.
"""
function read_constraints_from_yml(path::AbstractString)
    constraints_vec = Vector{Constraint}()
    open(path) do io
        yaml_iter = YAML.load_all(io, dicttype=OrderedDict)
        for d in yaml_iter
            pressure_GPa = d["pressure_kbar"] / 10
            temperature_C = d["temperature_C"]
            bulk = d["bulk"]
            bulk_oxides = d["bulk_oxides"]
            sys_in = d["sys_in"]
            assemblage = d["assemblage"]
            mineral_composition_apfu = d["mineral_comp_apfu"]
            mineral_elements = d["mineral_elements"]

            constraint = Constraint(pressure_GPa, temperature_C, bulk, bulk_oxides, sys_in, assemblage, mineral_composition_apfu, mineral_elements)
            push!(constraints_vec, constraint)
        end
    end
    return constraints_vec
end

"""
Generate constraints of the form P-T-Xbulk, that can be used when inverting a misfit to a functional relation (e.g., Ti-saturation in biotite)
In that case no mineral composition is needed, as the results of PEM are simply compared to a given functional relation.
"""
function gen_constraints_for_functional_inv(nb_constraints      ::Number;
                                            P_MIN_GPa           ::AbstractFloat                     = 0.2,
                                            P_MAX_GPa           ::AbstractFloat                     = 1.4,
                                            T_MIN_C             ::AbstractFloat                     = 450.,
                                            T_MAX_C             ::AbstractFloat                     = 700.,
                                            bulk_rock           ::Union{Vector{<:AbstractFloat}, Vector{Vector{<:AbstractFloat}}}   = [70.999, 0.758, 12.805, 6.342, 0.075, 3.978, 0.771, 1.481, 2.7895, 0.72933, 30],
                                            bulk_oxides         ::AbstractVector{<:AbstractString}  = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","O","H2O"],
                                            sys_in              ::AbstractString                    = "mol",
                                            phase               ::AbstractString                    = "bi",
                                            mineral_elements    ::AbstractVector{<:AbstractString}  = ["Si","Al","Ca","Mg","Fe","K","Na","Ti","O","Mn","H"],
                                            rand_seed           ::Union{Nothing, Integer}           = nothing)

    if !isnothing(rand_seed)
        rng = Xoshiro(rand_seed)
    else
        rng = Xoshiro()
    end

    P_GPa = rand(rng, nb_constraints) .* (P_MAX_GPa - P_MIN_GPa) .+ P_MIN_GPa
    T_C   = rand(rng, nb_constraints) .* (T_MAX_C - T_MIN_C)     .+ T_MIN_C

    # if only a single bulk vector is given, repeat this to use it for every constraint
    if typeof(bulk_rock) <: Vector{<:AbstractFloat}
        bulk_rock = repeat([bulk_rock], nb_constraints)
    end

    # generate the assemblage using MAGEMin_C for each constraint
    MAGEMin_db = Initialize_MAGEMin("mp", solver=2, verbose=false)

    out_vec = multi_point_minimization(P_GPa .* 10, T_C, MAGEMin_db, X=bulk_rock, Xoxides=bulk_oxides, sys_in=sys_in)
    assemblage = [out.ph for out in out_vec]

    Finalize_MAGEMin(MAGEMin_db)
    
    constraints_vec = Vector{Constraint}(undef, nb_constraints)

    for i in eachindex(1:nb_constraints)
        constraints_vec[i] = Constraint(P_GPa[i],
                                        T_C[i],
                                        bulk_rock[i],
                                        bulk_oxides,
                                        sys_in,
                                        assemblage[i],
                                        nothing,
                                        mineral_elements)
    end

    return constraints_vec
end


"""
to do..
"""
function read_bulks()
    #//TODO - read fucntions to import bulks form the FPWMP database (read from a csv)
    return bulk_vec
end
