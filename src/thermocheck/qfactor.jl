#=
Assess if a set thermodynamic data reproduces the mineral assemblage (Qasm) and modes (Qvol) for selected P-T points
and bulk rock composition (corresponding to "key pseudosections/MADs") with respect to the natural data or in to other
thermodynamic databases.
=#
using ..ThermoFit

const Pkbar_TC_POINTS = [
    (2.0, 550.),
    (2.0, 600.),
    (2.0, 650.),
    (2.0, 700.),
    (5.0, 550.),
    (5.0, 600.),
    (5.0, 650.),
    (5.0, 700.),
    (10.0, 550.),
    (10.0, 600.),
    (10.0, 650.),
    (10.0, 700.)
]

const BULK_OXIDES = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "O", "H2O"]

const BULK_ROCKS = (
    FPWMP         = [70.999, 0.758, 12.805, 6.342, 0.075, 3.978, 0.771, 1.481, 2.7895, 0.72933, 30.0],
    # FPWMP_projCa  = [],
    # Fig7_W14b     = []
)

"""
Reference assemblages from [Pattison & Forshaw (2025)](https://doi.org/10.1093/petrology/egaf039) Fig 2.
All assemblages contain Ti-Oxides (Ilm/Rt) on top of the listed ones.
"""
const ASSEMBLAGE_PF25_PETROGENETIC_GRID = [
    ["mu", "cd", "bi", "pl", "q"],                  # 2.0 kbar - 550 °C
    ["mu", "cd", "bi", "pl", "q"],                  # 2.0 kbar - 600 °C
    ["afs", "cd", "bi", "pl", "q"],                 # 2.0 kbar - 650 °C
    ["afs", "cd", "bi", "pl", "q", "liq"],          # 2.0 kbar - 700 °C
    ["mu", "chl", "g", "bi", "pl", "q"],            # 5.0 kbar - 550 °C
    ["st", "±sil", "mu", "g", "bi", "pl", "q"],     # 5.0 kbar - 600 °C
    ["sil", "mu", "g", "bi", "pl", "q"],            # 5.0 kbar - 650 °C
    ["sil", "g", "afs", "liq", "bi", "pl", "q"],    # 5.0 kbar - 700 °C
    ["mu", "chl", "g", "bi", "pl", "q"],            # 10.0 kbar - 550 °C
    ["mu", "g", "bi", "pl", "q"],                   # 10.0 kbar - 600 °C
    ["mu", "g", "ky", "bi", "pl", "q"],             # 10.0 kbar - 650 °C
    ["mu", "g", "ky", "liq", "bi", "pl", "q"]       # 10.0 kbar - 700 °C
]


"""
Q_asm (Quality factor for the mineral assemblage) after [Duesterhoeft & Lanari (2020)](https://doi.org/10.1111/jmg.12538)
"""
function q_asm(asm_observed  ::Vector{<:AbstractString},
               asm_reference ::Vector{<:AbstractString})
    
    l = length(intersect(asm_observed, asm_reference))
    k = length(union(asm_observed, asm_reference))
    
    return 100 * (l/k)
end


"""
Q_vol (Quality factor for the mineral modes) after [Duesterhoeft & Lanari (2020)](https://doi.org/10.1111/jmg.12538)
"""
function q_vol(asm_observed  ::Vector{<:AbstractString},
               asm_reference ::Vector{<:AbstractString},
               modes_observed::Vector{<:Number},
               modes_reference::Vector{<:Number})
    
    l = intersect(asm_observed, asm_reference)

    # filter and reorder the modes to match the order and reduce to the "l" matching phases
    filter_obs = filter!(!isnothing, [findfirst(.==(phase), l) for phase in asm_observed])
    filter_ref = filter!(!isnothing, [findfirst(.==(phase), l) for phase in asm_reference])
    modes_observed = modes_observed[filter_obs]
    modes_reference = modes_reference[filter_ref]

    return 100 * sqrt(sum((modes_observed .+ modes_reference) ./ 2 .* (1 .- abs.(modes_observed .- modes_reference) ./ (modes_observed .+ modes_reference)).^2)) 
end


"""
Create test points (<:Constraint) for each P-T point in `Pkbar_TC_POINTS` for the specified bulk rock composition
"""
function create_test_points(bulk_name ::Union{AbstractString, Symbol})
    # ensure the bulk_name is a Symbol
    if isa(bulk_name, AbstractString)
        bulk_name = Symbol(bulk_name)
    end
    bulk = BULK_ROCKS[bulk_name]

    constraints = Vector{Constraint}(undef, length(Pkbar_TC_POINTS))
    
    for i in eachindex(Pkbar_TC_POINTS)
        p, t = Pkbar_TC_POINTS[i]
        
        constraints[i] = Constraint(p * 0.1, 
                                    t,
                                    bulk,
                                    BULK_OXIDES,
                                    "mol",
                                    nothing,
                                    nothing,
                                    ["Si", "Al", "Ca", "Mg", "Fe", "K", "Na", "Ti", "O", "Mn", "H"])
    end
    return constraints
end


"""
Call the forward model (G-minimisation) for a set of thermodynamic (database with optionally modified Ws and G°-corrections)
to probe the predictions of the mineral assemblage and modes for a set of P-T-Xbulk test points (<:Constraint).
"""
function probe_test_points(test_points ::Vector{<:ThermoFit.Constraint},
                           database    ::AbstractString,
                           phase       ::AbstractString;
                           w_g         ::Union{Nothing, Vector{<:Number}} = nothing,
                           g0_corr_vec ::Union{Nothing, Vector{<:Number}} = nothing,
                           g0_corr_em  ::Union{Nothing, Vector{<:AbstractString}} = nothing)

    MAGEMin_db = Initialize_MAGEMin(database, solver=2, verbose=false)

    asm_vec = Vector{Vector{<:AbstractString}}(undef, length(test_points))
    modes_vec = Vector{Vector{<:Number}}(undef, length(test_points))

    @threads for i in eachindex(test_points)
        constraint = test_points[i]
            
        # identify thread and acess the MAGEMin_db of the thread
        id          = Threads.threadid()

        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        out = forward_call(phase, database, constraint, gv, z_b, DB, splx_data, w_g = w_g, g0_corr_vec = g0_corr_vec, g0_corr_em = g0_corr_em)

        # extract the mineral assemblage and modes from the output
        asm_vec[i] = copy(out.ph)
        modes_vec[i] = copy(out.ph_frac_vol)
    end

    Finalize_MAGEMin(MAGEMin_db)
    return asm_vec, modes_vec
end


function plot_qasm(database       ::AbstractString,
                   reference      ::String,
                   phase          ::AbstractString;
                   w_g            ::Union{Nothing, Vector{<:Number}} = nothing,
                   g0_corr_vec    ::Union{Nothing, Vector{<:Number}} = nothing,
                   g0_corr_em     ::Union{Nothing, Vector{<:AbstractString}} = nothing)

    if reference == "PF25"
        asm_vec_ref = ASSEMBLAGE_PF25_PETROGENETIC_GRID
    elseif reference == "W14"
        nothing
    else
        error("Unknown reference: $reference")
    end

    fig = Figure(size=(600, length(BULK_ROCKS)*400))
    p = [p for (p, _) in Pkbar_TC_POINTS]
    t = [t for (_, t) in Pkbar_TC_POINTS]

    for (i, bulk_name) in enumerate(keys(BULK_ROCKS))
        test_points = create_test_points(bulk_name)
        asm_vec, _ = probe_test_points(test_points, database, phase; w_g = w_g, g0_corr_vec = g0_corr_vec, g0_corr_em = g0_corr_em)
        println(asm_vec)
        if reference == "PF25"
            # strip ilm and/or rt from the test points' assemblages
            # these are considered as Ti-oxides that can always be present
            # in the ref assemblage of PF25
            asm_vec = [filter(x -> x != "ilm" && x != "ru", asm) for asm in asm_vec]
            # strip the H20 as this is not considered in the PF25 reference assemblages
            asm_vec = [filter(x -> x != "H2O", asm) for asm in asm_vec]
        elseif reference == "W14"
            ref_points = create_test_points(bulk_name)
            asm_vec_ref, _ = probe_test_points(ref_points, "mp", phase)
        end
        
        q_asm_vec = Vector{Float64}(undef, length(asm_vec))
        for i in eachindex(asm_vec)
            q_asm_vec[i] = q_asm(asm_vec[i], asm_vec_ref[i])
        end
        # calculate the differences between the observed and reference assemblages
        asm_diff_vec = asm_diff_vector(asm_vec, asm_vec_ref)

        ax = Axis(fig[i, 1], title = "Bulk: $bulk_name")
        s = scatter!(ax, t, p, color = q_asm_vec, colormap = :batlow, colorrange = (0, 100), markersize = 60)
        text!(ax, t, p, text=string.(round.(Int, q_asm_vec)),
              align=(:center, :center), font=:bold)
        text!(ax, t, p.-0.85, text=asm_diff_vec,
              align=(:center, :bottom), fontsize=12)

        Colorbar(fig[i, 2], s, label = "Q_asm(model - reference) [%]", width = 20, height = 300)
        ax.xlabel = "Temperature [°C]"
        ax.ylabel = "Pressure [kbar]"
        xlims!(ax, 500, 750)
        ylims!(ax, 0, 12)
    end
    return fig
end


function asm_diff_vector(asm_vec,
                         asm_ref)
    
    asm_diff_vec = Vector{String}(undef, length(asm_vec))

    for (i, (asm, asm_ref)) in enumerate(zip(asm_vec, asm_ref))
        # find all strings present in asm but not in asm_ref
        phase_added = setdiff(asm, asm_ref)
        phases_added = ["+" * phase for phase in phase_added]
        # find all strings present in asm_ref but not in asm
        phase_missing = setdiff(asm_ref, asm)
        phases_missing = ["-" * phase for phase in phase_missing]
        # combine the added and missing phases
        asm_diff = vcat(phases_added, phases_missing)
        # sort the phases alphabetically
        sort!(asm_diff)
        # convert to a string
        asm_diff = join(asm_diff, ", ")
        asm_diff_vec[i] = asm_diff == "" ? "no diff" : asm_diff
    end

    return asm_diff_vec
end
