struct Constraint
    pressure_GPa::Float64
    temperature_C::Float64
    bulk::Vector{Float64}
    bulk_oxides::Vector{String}
    sys_in::String
    mineral_composition_apfu::Vector{Float64}
    mineral_elements::Vector{String}
end


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
            mineral_compositions_apfu = d["mineral_compositions_apfu"]
            println(mineral_compositions_apfu)

            constraint = Constraint(pressure_GPa, temperature_C, bulk, bulk_oxides, sys_in, mineral_compositions_apfu)
            push!(constraints_vec, constraint)
        end
    end
    return constraints_vec
end