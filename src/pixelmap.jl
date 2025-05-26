
"""
    pixelmap(temperature_vec, pressure_vec, bulk, bulk_oxides, phase, database, w_g, sys_in)
Pixelmap-like caluclation using MAGEMin for quick visualisation of the effect of altered W_G and G0.
"""
function pixelmap(temperature_vec          ::Vector{Float64},
                  pressure_vec             ::Vector{Float64},
                  bulk                     ::Vector{Float64},
                  bulk_oxides              ::AbstractVector{String},
                  database                 ::String,
                  sys_in                   ::String,
                  comp_variables_export    ::AbstractVector{String},
                  phase                    ::String;
                  w_g                      ::Union{Vector{Float64}, Nothing} = nothing,
                  G_0                      ::Union{Vector{Float64}, Nothing} = nothing)

    # create and populate matrices for P and T
    pressure_matrix    = zeros(Float64, length(pressure_vec), length(temperature_vec))
    temperature_matrix = zeros(Float64, length(pressure_vec), length(temperature_vec))
    
    for i in axes(pressure_matrix, 1), j in axes(pressure_matrix, 2)
        pressure_matrix[i, j] = pressure_vec[i]
        temperature_matrix[i, j] = temperature_vec[j]
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
