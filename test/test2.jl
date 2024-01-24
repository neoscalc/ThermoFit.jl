using MAGEMin_C
using NLopt

database = "mp";

global gv, z_b, DB, splx_data  = init_MAGEMin(database);


P          = 5
T          = 650
gv = use_predefined_bulk_rock(gv, 0, database);
# initialize MAGEMin up to G0 and W's point
gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);


# get the solution phase structure (size gv.len_ss)
ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

print("\n---------------------------------\n");
print("    Database information (",database,")\n");
print("----------------------------------\n");
for i=1:gv.len_ss
    print("   ",ss_names[i],": ",i,"; n_W's, ",ss_struct[i].n_w,"; n_em's, ",ss_struct[i].n_em,"\n")
end
print("----------------------------------\n\n");

W = unsafe_wrap(Vector{Cdouble},ss_struct[3].W, ss_struct[3].n_w)
W[2] = 24.0

out       = pwm_run(gv, z_b, DB, splx_data);