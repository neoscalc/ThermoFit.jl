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
function read_constraints_from_yml(path::AbstractString; log_io::IO = stdout)
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

    print_constraints(length(constraints_vec),
                      constraints_yaml=path,
                      io=log_io)
    return constraints_vec
end


function save_constraints_to_yml(constraints_vec::Vector{Constraint}, path::AbstractString; DATABASE="mp")

    doc_keywords = ["GEN DATA $(i)" for i in 1:length(constraints_vec)]

    open(path, "w") do io
        for (i, constraint) in enumerate(constraints_vec)
            write(io, "---\n")
            write(io, "# $(doc_keywords[i])\n")
            YAML.write(io,
                       OrderedDict(
                            "pressure_kbar" => constraint.pressure_GPa * 10,
                            "temperature_C" => constraint.temperature_C,
                            "bulk" => constraint.bulk,
                            "bulk_oxides" => constraint.bulk_oxides,
                            "sys_in" => constraint.sys_in,
                            "database" => DATABASE,
                            "assemblage" => constraint.assemblage,
                            "mineral_comp_apfu" => constraint.mineral_composition_apfu,
                            "mineral_elements" => constraint.mineral_elements))
        end
    end
    
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
                                            bulk_rock           ::Union{AbstractVector, AbstractString} = [70.999, 0.758, 12.805, 6.342, 0.075, 3.978, 0.771, 1.481, 2.7895, 0.72933, 30],
                                            bulk_oxides         ::AbstractVector{<:AbstractString}  = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","O","H2O"],
                                            sys_in              ::AbstractString                    = "mol",
                                            phase               ::AbstractString                    = "bi",
                                            mineral_elements    ::AbstractVector{<:AbstractString}  = ["Si","Al","Ca","Mg","Fe","K","Na","Ti","O","Mn","H"],
                                            rng                 ::Union{Union{Nothing, Integer}, AbstractRNG} = nothing,
                                            save_to_yaml        ::Bool                              = false,
                                            log_io              ::IO = stdout)

    if isnothing(rng)
        rng = Xoshiro()
    elseif typeof(rng) <: Integer
        rng = Xoshiro(rng)
    end

    if typeof(bulk_rock) <: AbstractString
        bulk_rocks, bulk_oxides = read_FPWMP_bulks(bulk_rock)
        sys_in = "wt%"
    end

    MAGEMin_db = Initialize_MAGEMin("mp", solver=2, verbose=false)
    constraints_vec = Vector{Constraint}(undef, nb_constraints)
    nb_constraints_generated = 0
    while nb_constraints_generated < nb_constraints
        P_GPa = rand(rng) .* (P_MAX_GPa - P_MIN_GPa) .+ P_MIN_GPa
        T_C   = rand(rng) .* (T_MAX_C - T_MIN_C)     .+ T_MIN_C

        # if a single bulk vector is given use this open
        # if not, random sample from the loaded FPWMP_bulks 
        if typeof(bulk_rock) <: AbstractString
            bulk_rock = rand(rng, bulk_rocks)
        end

        out = single_point_minimization(P_GPa .* 10, T_C, MAGEMin_db, X=bulk_rock, Xoxides=bulk_oxides, sys_in=sys_in)
        assemblage = out.ph

        if phase in assemblage
            comp_apfu = out.SS_vec[findfirst(phase .== assemblage)].Comp_apfu
            mineral_composition_apfu = OrderedDict{String, Any}()
            mineral_composition_apfu[phase] = comp_apfu

            i = nb_constraints_generated + 1
            constraints_vec[i] = Constraint(P_GPa,
                                            T_C,
                                            bulk_rock,
                                            bulk_oxides,
                                            sys_in,
                                            assemblage,
                                            mineral_composition_apfu,
                                            mineral_elements)
            
            nb_constraints_generated += 1
        end
    end

    Finalize_MAGEMin(MAGEMin_db)

    print_constraints(nb_constraints,
                      constraints_yaml = "generated",
                      constraints_gen  = true,
                      P_MIN_GPa        = P_MIN_GPa,
                      P_MAX_GPa        = P_MAX_GPa,
                      T_MIN_C          = T_MIN_C,
                      T_MAX_C          = T_MAX_C,
                      bulk_rocks        = string(bulk_rock),
                      sys_in           = sys_in,
                      io               = log_io)

    if save_to_yaml
        date = Dates.format(Dates.now(), "yyyy-mm-dd_HHMM")
        save_constraints_to_yml(constraints_vec, "generated_constraints_$date.yml")
    end
    return constraints_vec
end


"""
Read bulk rock composition (in wt%) from a csv of the Forshaw & Pattison worldwide metapelitic whole rock database.

Optionally, the bulk rock can be renormalised to 100 wt% and/or the P2O5 can be projected from apatite to correct the CaO content.
"""
function read_FPWMP_bulks(file_path             ::AbstractString;
                          buffer_H2O            ::Float64 = 30.0,
                          project_from_apatite  ::Bool = false,
                          renormalise_anhydrous ::Bool = false)

    bulks_df = CSV.read(file_path, DataFrame)

    if project_from_apatite
        bulks_df = DataFrame(map(ThermoFit.project_from_apatite, eachrow(bulks_df)))

    else # drop P2O5
        bulks_df = select(bulks_df, Not("P2O5"))
    end

    # replace missing values with 0
    bulks_df = coalesce.(bulks_df, 0.0)

    bulk = [collect(b) for b in eachrow(bulks_df)]
    
    if renormalise_anhydrous
        # renormalise the anhydrous bulk to 100 wt%
        bulk = [b ./ sum(b) .* 100 for b in bulk]
    end

    # add H2O to the bulk
    for i in eachindex(bulk)
        bulk[i] = vcat(bulk[i], buffer_H2O)
    end

    bulk_oxides = names(bulks_df)
    bulk_oxides =vcat(bulk_oxides, "H2O")
        
    return bulk, bulk_oxides
end

function project_from_apatite(bulk  ::DataFrameRow)
    MOLARMASS_P2O5 = 283.88
    MOLARMASS_CaO  = 56.08
    
    CaO_apatite = 10/3 * bulk["P2O5"] * MOLARMASS_CaO / MOLARMASS_P2O5

    # remove the CaO_apatite from the CaO
    if ismissing(CaO_apatite)
        CaO_apatite = 0    
    elseif bulk["CaO"] > CaO_apatite
        bulk["CaO"] -= CaO_apatite
    else 
        @warn "CaO_apatite is larger than the bulk CaO, setting it to 0"
        bulk["CaO"] = 0
    end

    # remove the P2O5 from the bulk
    bulk = bulk[Not("P2O5")]
    return bulk
end
