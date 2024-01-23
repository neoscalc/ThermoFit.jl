using MAGEMin_C

# Initialize database  - new way
db          =   "mp"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b); um, ultramafic (Evans & Frost 2021)
data        =   Initialize_MAGEMin(db, verbose=true)

P           =   4.0
T           =   500.0

out     =   single_point_minimization(P, T, data, test=0)

out.SS_vec[6]