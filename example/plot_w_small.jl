# Test of the function
using ThermoFit
# using Base.Threads

# define the bulk composition
bulk = [70.999, 12.805, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0]
bulk_oxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"]

# define steps between 450 and 800 °C
resolution = 30
ti = Vector(range(400, stop=800, length=resolution))
pi = Vector(range(1, stop=12, length=resolution))

# define other parameters
database = "mp"
phase = "bi"
sys_in = "mol"
w_g = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]

# define the variables
variables_export = ["Si", "Ti", "Fe", "Mg", "Mn"]

# call the function
result = pixelmap_calc(ti, pi, bulk, bulk_oxides, phase, database, w_g, sys_in, variables_export)

# w_g_new = [12.0, 4.0, 10.0, 30.0, 8.0, 9.0, 8.0, 15.0, 32.0, 13.6, 6.3, 7.0, 24.0, 5.6, 8.1, 40.0, 1.0, 13.0, 40.0, 30.0, 11.6]
#w_g_new = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

w_g_new = [9.3, 10.6, 17.8, 19.15, 15.24, 13.64, 10.50, 14.17, 38.32, 17.18, 19.57, 8.35, 1.04, 11.13, 16.83, 22.97, 12.99, 16.34, 19.43, 14.59, 14.40]

# call the function
result_new = pixelmap_calc(ti, pi, bulk, bulk_oxides, phase, database, w_g_new, sys_in, variables_export)

# plot the results from columns variables_export using heatmaps (one subplot per variable)
using Plots

si_lim = (2.4,3)
ti_lim = (0,0.2)
fe_lim = (0.5,2)
mg_lim = (0.5,2)
mn_lim = (0,0.03)
xmg_lim = (0.3,0.8)

Si = reshape(result.Si, resolution, resolution)
p1 = heatmap(ti, pi, Si, title="Si", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = si_lim)

Si_new = reshape(result_new.Si, resolution, resolution)
p2 = heatmap(ti, pi, Si_new, title="Si_new", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = si_lim)

Ti = reshape(result.Ti, resolution, resolution)
p3 = heatmap(ti, pi, Ti, title="Ti", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = ti_lim)

Ti_new = reshape(result_new.Ti, resolution, resolution)
p4 = heatmap(ti, pi, Ti_new, title="Ti_new", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = ti_lim)

Fe = reshape(result.Fe, resolution, resolution)
p5 = heatmap(ti, pi, Fe, title="Fe", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = fe_lim)

Fe_new = reshape(result_new.Fe, resolution, resolution)
p6 = heatmap(ti, pi, Fe_new, title="Fe_new", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = fe_lim)

Mg = reshape(result.Mg, resolution, resolution)
p7 = heatmap(ti, pi, Mg, title="Mg", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = mg_lim)

Mg_new = reshape(result_new.Mg, resolution, resolution)
p8 = heatmap(ti, pi, Mg_new, title="Mg_new", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = mg_lim)

Mn = reshape(result.Mn, resolution, resolution)
p9 = heatmap(ti, pi, Mn, title="Mn", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = mn_lim)

Mn_new = reshape(result_new.Mn, resolution, resolution)
p10 = heatmap(ti, pi, Mn_new, title="Mn_new", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = mn_lim)

XMg = Mg ./ (Mg .+ Fe)
p11 = heatmap(ti, pi, XMg, title="XMg", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = xmg_lim)

XMg_new = Mg_new ./ (Mg_new .+ Fe_new)
p12 = heatmap(ti, pi, XMg_new, title="XMg_new", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = xmg_lim)

plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, layout=(6,2), size=(1000, 1500))



# Plot differences in Si
heatmap(ti, pi, Si_new-Si, title="Si_new-Si", xlabel="T (°C)", ylabel="P (kbar)", c=:viridis, clim = (-0.05,0.05))

