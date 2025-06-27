
"""
    pixelmap(temperature_vec, pressure_vec, bulk, bulk_oxides, phase, database, w_g, sys_in)
Pixelmap-like caluclation using MAGEMin for quick visualisation of the effect of altered W_G and G0.
"""
function pixelmap(temperature_C_vec        ::Vector{Float64},
                  pressure_kbar_vec        ::Vector{Float64},
                  bulk                     ::Vector{Float64},
                  bulk_oxides              ::AbstractVector{String},
                  database                 ::String,
                  sys_in                   ::String,
                  comp_variables_export    ::AbstractVector{String},
                  phase                    ::String;
                  w_g                      ::Union{Vector{Float64}, Nothing} = nothing,
                  G_0                      ::Union{Vector{Float64}, Nothing} = nothing)

    # create and populate matrices for P and T
    pressure_matrix    = zeros(Float64, length(pressure_kbar_vec), length(temperature_C_vec))
    temperature_matrix = zeros(Float64, length(pressure_kbar_vec), length(temperature_C_vec))
    
    for i in axes(pressure_matrix, 1), j in axes(pressure_matrix, 2)
        pressure_matrix[i, j] = pressure_kbar_vec[i]
        temperature_matrix[i, j] = temperature_C_vec[j]
    end
    
    # create vectors (P, T, mineral comp) for the calculation 
    p_vec = vec(pressure_matrix)
    t_vec = vec(temperature_matrix)
    min_comp = Array{Union{Float64, Missing}}(missing, length(p_vec), length(comp_variables_export))
    
    MAGEMin_db  = Initialize_MAGEMin(database, verbose=false, solver=2)
    
    @threads for i in ProgressBar(eachindex(p_vec))
        # identify thread and acess the MAGEMin_db of the thread
        id          = Threads.threadid()
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]
    
        # init minimisation
        gv = define_bulk_rock(gv, bulk, bulk_oxides, sys_in, database)
        gv, z_b, DB, splx_data = pwm_init(p_vec[i], t_vec[i], gv, z_b, DB, splx_data)
    
        # alter W_G and G0
        ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss))
        ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss)
    
        ss_idx = findfirst(ss_names .== phase)
        if !isnothing(w_g)
            unsafe_copyto!(ss_struct[ss_idx].W, pointer(w_g), ss_struct[ss_idx].n_w)
        end
        if !isnothing(G_0)
            g_base = unsafe_wrap(Vector{Float64}, ss_struct[ss_idx].gbase, ss_struct[ss_idx].n_em)
            g_base .+= G_0
            unsafe_copyto!(ss_struct[ss_idx].gbase, pointer(g_base), ss_struct[ss_idx].n_em)
        end
    
        # run minimisation
        out = pwm_run(gv, z_b, DB, splx_data)
    
        # extract mineral composition, only the compositional varibales specified in comp_variables_export
        if phase in out.ph
            idx_elem_export_in_out = indexin(comp_variables_export, out.elements)
        
            min_comp_i = out.SS_vec[findfirst(out.ph .== phase)].Comp_apfu
            min_comp_i = min_comp_i[idx_elem_export_in_out]
        
            min_comp[i, :] .= min_comp_i
        end
    end
    
    Finalize_MAGEMin(MAGEMin_db)
    
    # shape the mineral composition matrix to the same shape as the pressure_matrix/temperature_matrix
    # mineral_composition is a vector of matrices, each matrix holds the compositional variable values over P and T
    mineral_composition = []
    for i in eachindex(comp_variables_export)
        min_comp_matrix = reshape(min_comp[:, i], size(pressure_matrix))
        push!(mineral_composition, min_comp_matrix)
    end
    
    return mineral_composition, pressure_matrix, temperature_matrix
end


"""
Plot figure with pixel maps of Ti-in-Bt temperatures and ϕ(Ti-in-Bt - Tmodel) misfit comparing original and modified thermodynamic parameters.
"""
function pixmap_misfit_Ti_in_Bt_temperatures(file_path                ::String;
                                             w_g                      ::Union{Vector{Float64}, Nothing} = nothing,
                                             G_0                      ::Union{Vector{Float64}, Nothing} = nothing,
                                             temperature_C_vec        ::Vector{Float64}         = Vector(450.:25:750.),
                                             pressure_kbar_vec        ::Vector{Float64}         = Vector(2.:0.5:14.),
                                             bulk                     ::Vector{Float64}         = [70.999, 0.758, 12.805, 6.342, 0.075, 3.978, 0.771, 1.481, 2.7895, 0.72933, 30.0],
                                             bulk_oxides              ::AbstractVector{String}  = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "O", "H2O"],
                                             database                 ::String                  = "mp",
                                             sys_in                   ::String                  = "mol",
                                             comp_variables_export    ::AbstractVector{String}  = ["Si", "Al", "Ca", "Mg", "Fe", "K", "Na", "Ti", "O", "Mn", "H"],
                                             phase                    ::String                  = "bi",
                                             colormap_Tmap            ::String                  = "berlin",
                                             colormap_misfitmap       ::String                  = "bam")

    # pixel map with the modified W_G and G0
    pxm, p_, t_ = pixelmap(temperature_C_vec,
                           pressure_kbar_vec,
                           bulk,
                           bulk_oxides,
                           database,
                           sys_in,
                           comp_variables_export,
                           phase,
                           w_g = w_g,
                           G_0 = G_0)

    # pixel map with the original W_G and G0
    pxo, p_, t_ = pixelmap(temperature_C_vec,
                           pressure_kbar_vec,
                           bulk,
                           bulk_oxides,
                           database,
                           sys_in,
                           comp_variables_export,
                           phase)

    # predict temperature using Ti-in-Bt thermometry (Henry et al. 2005)
    pxm_TBt05 = ThermoFit.Ti_in_Bt_Henry05.(pxm[findfirst("Ti" .== comp_variables_export)], pxm[findfirst("Mg" .== comp_variables_export)], pxm[findfirst("Fe" .== comp_variables_export)])
    pxo_TBt05 = ThermoFit.Ti_in_Bt_Henry05.(pxo[findfirst("Ti" .== comp_variables_export)], pxo[findfirst("Mg" .== comp_variables_export)], pxo[findfirst("Fe" .== comp_variables_export)])
    # calculate msifit between Ti-in-Bt temperature and the temperature used for PEM
    misfit_o = pxo_TBt05 .- t_
    misfit_m = pxm_TBt05 .- t_

    # extrema are calculated for the colorrange in the plotting
    extrema_pxo_TBt05 = extrema(skipmissing(pxo_TBt05))
    extrema_pxm_TBt05 = extrema(skipmissing(pxm_TBt05))

    clims = (min(first(extrema_pxo_TBt05), first(extrema_pxm_TBt05)),
             max(last(extrema_pxo_TBt05), last(extrema_pxo_TBt05)))

    # center the clims_misfit so that they are symmetric around zero otherwise,
    # the relatively small/rare overestimation of temps by Henry's are visually over represented
    extrema_misfit = extrema(skipmissing(vcat(misfit_o, misfit_m)))

    clims_misfit = (-max(abs(first(extrema_misfit)), abs(last(extrema_misfit))),
                    max(abs(first(extrema_misfit)), abs(last(extrema_misfit))))

    fig = Figure(size=(800, 800))

    ax1 = Axis(fig[1, 1])
    ax1.title = "T [°C] Ti-in-Bt - original Ws"
    ax1.aspect = 1.0

    hm = heatmap!(ax1, temperature_C_vec, pressure_kbar_vec, pxo_TBt05', colormap=colormap_Tmap, colorrange=clims)
    Colorbar(fig[1, 2], colormap=colormap_Tmap, colorrange=clims)

    ax1 = Axis(fig[1, 3])
    ax1.title = "ΔT(Ti-in-Bt - Tmodel) - original Ws"
    ax1.aspect = 1.0

    hm = heatmap!(ax1, temperature_C_vec, pressure_kbar_vec, misfit_o', colormap=colormap_misfitmap, colorrange=clims_misfit)
    Colorbar(fig[1, 4], colormap=colormap_misfitmap, colorrange=clims_misfit)

    ax2 = Axis(fig[2, 1])
    ax2.title = "T [°C] Ti-in-Bt - modified Ws"
    ax2.aspect = 1.0

    hm = heatmap!(ax2, temperature_C_vec, pressure_kbar_vec, pxm_TBt05', colormap=colormap_Tmap, colorrange=clims)
    Colorbar(fig[2, 2], colormap=colormap_Tmap, colorrange=clims)

    ax2 = Axis(fig[2, 3])
    ax2.title = "ΔT(Ti-in-Bt - Tmodel) - modified Ws"
    ax2.aspect = 1.0

    hm = heatmap!(ax2, temperature_C_vec, pressure_kbar_vec, misfit_m', colormap=colormap_misfitmap, colorrange=clims_misfit)
    Colorbar(fig[2, 4], colormap=colormap_misfitmap, colorrange=clims_misfit)

    save(file_path, fig)

    return fig
end


"""
Plot pixel maps of each element in <PHASE> comparing original and modified thermodynamic parameters.
"""
function pixmap_phase_comp_og_v_mod(file_path                ::String;
                                    w_g                      ::Union{Vector{Float64}, Nothing} = nothing,
                                    G_0                      ::Union{Vector{Float64}, Nothing} = nothing,
                                    temperature_C_vec        ::Vector{Float64}         = Vector(450.:25:750.),
                                    pressure_kbar_vec        ::Vector{Float64}         = Vector(2.:0.5:14.),
                                    bulk                     ::Vector{Float64}         = [70.999, 0.758, 12.805, 6.342, 0.075, 3.978, 0.771, 1.481, 2.7895, 0.72933, 30.0],
                                    bulk_oxides              ::AbstractVector{String}  = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "O", "H2O"],
                                    database                 ::String                  = "mp",
                                    sys_in                   ::String                  = "mol",
                                    comp_variables_export    ::AbstractVector{String}  = ["Si", "Al", "Ca", "Mg", "Fe", "K", "Na", "Ti", "O", "Mn", "H"],
                                    phase                    ::String                  = "bi",
                                    colormap                 ::String                  = "batlow")


    # pixel map with the modified W_G and G0
    pxm, p_, t_ = pixelmap(temperature_C_vec,
                           pressure_kbar_vec,
                           bulk,
                           bulk_oxides,
                           database,
                           sys_in,
                           comp_variables_export,
                           phase,
                           w_g = w_g,
                           G_0 = G_0)

    # pixel map with the original W_G and G0
    pxo, p_, t_ = pixelmap(temperature_C_vec,
                           pressure_kbar_vec,
                           bulk,
                           bulk_oxides,
                           database,
                           sys_in,
                           comp_variables_export,
                           phase)
    # extrema are calculated for the colorrange in the plotting
    extrema_pxo = [extrema(skipmissing(x)) for x in pxo]
    extrema_pxm = [extrema(skipmissing(x)) for x in pxm]

    fig = Figure(size=(800, length(comp_variables_export)*400))
    for i in eachindex(comp_variables_export)
         ax = Axis(fig[i, 1])

         # define clims
         clims = (min(first(extrema_pxo[i]), first(extrema_pxm[i])),
                  max(last(extrema_pxo[i]), last(extrema_pxm[i]))+ 0.001)  # add small offset to avoid (0,0)


         e = comp_variables_export[i]
         ax.title = "$e in Bi [apfu] - original Ws"
         ax.aspect = 1.0
         hm = heatmap!(ax, temperature_C_vec, pressure_kbar_vec, pxo[i]', colormap=colormap, colorrange=clims)

         ax = Axis(fig[i, 2])

         e = comp_variables_export[i]
         ax.title = "$e in Bi [apfu] - modified Ws"
         ax.aspect = 1.0
         hm = heatmap!(ax, temperature_C_vec, pressure_kbar_vec, pxm[i]', colormap=colormap, colorrange=clims)

         Colorbar(fig[i, 3], colormap=colormap, colorrange=clims)
    end
    save(file_path, fig)

    return fig
end
