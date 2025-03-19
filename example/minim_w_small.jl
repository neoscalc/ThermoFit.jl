# ---------------------------------------------------------------------------------------------
#                    ------------------------------------------------
#                   |         THERMOFIT 0.0.1 - EXAMPLE FILE         |
#                    ------------------------------------------------
#
# This file shows an example of how to use the ThermoFit package to adjust Margules parameters.
# 
# As a proof of concept the example presented here demonstrates how the thermodynamic parameters of the biotite solid solution 
# defined in White et al. (2014) can be retireved from altered parameeters using inversion modelling with ThermoFit.
# As constraints 49 P-T-X + biotite composition data points simulated with MAGEMin using the White et al. 2014 solution model are used.
#
# Philip Hartmeier & Pierre Lanari
# February 2024 (Cassis, France)
#
#        | /=   `/ () |_| /?   ~|~ |-| [- /? |\/| () |) `/ |\| /\ |\/| | (   |) /\ ~|~ /\
#                 _\~ |_| ( /< _\~   |_| _\~ [-   ~|~ |-| [- /? |\/| () /= | ~|~
#
# ---------------------------------------------------------------------------------------------

# 1. Change the working directory to the folder where this file is located
cd(@__DIR__)

# 2. Import ThermoFit.jl
using ThermoFit

# 3. Initialise the data for the inversion
thermodynamic_database = "mp";
solid_solution = "bi";

# start with the Margules for biotite from White et al. (2014)
# alter the Margules W(phl,annm) 12 -> 6 kJ/mol
# alter the Margules W(phl,obi)   4 -> 2  kJ/mol
# alter the Margules W(phl,east) 10 -> 20 kJ/mol
# allow the altered Margules to change within a range of 0 - 50 kJ/mol
w_names =  ["W(phl,annm)",
            "W(phl,obi)",
            "W(phl,east)",
            "W(phl,tbi)",
            "W(phl,fbi)",
            "W(phl,mmbi)",
            "W(annm,obi)",
            "W(annm,east)",
            "W(annm,tbi)",
            "W(annm,fbi)",
            "W(annm,mmbi)",
            "W(obi,east)",
            "W(obi,tbi)",
            "W(obi,fbi)",
            "W(obi,mmbi)",
            "W(east,tbi)",
            "W(east,fbi)",
            "W(east,mmbi)",
            "W(tbi,fbi)",
            "W(tbi,mmbi)",
            "W(fbi,mmbi)"];

w_initial_values = [6    0  0 ;
                    2    0  0 ;
                    20   0  0 ;
                    30   0  0 ;
                    8    0  0 ;
                    9    0  0 ;
                    8    0  0 ;
                    15   0  0 ;
                    32   0  0 ;
                    13.6 0  0 ;
                    6.3  0  0 ;
                    7    0  0 ;
                    24   0  0 ;
                    5.6  0  0 ;
                    8.1  0  0 ;
                    40   0  0 ;
                    1    0  0 ;
                    13   0  0 ;
                    40   0  0 ;
                    30   0  0 ;
                    11.6 0  0];

w_lower_bounds = copy(w_initial_values);
w_lower_bounds[1:3,1] .= 0.;
w_upper_bounds = copy(w_initial_values);
w_upper_bounds[1:3,1] .= 50.;

# use a G0 correction for the endmember annm of biotite
# Increase the G0 of annm by 2 kJ/mol, reversing the ann -> annm G0 correction
g0_corr_endmembers      = ["annm"];
g0_corr_initial_values  = [2.0];
g0_corr_lower_bounds    = [-3.];
g0_corr_upper_bounds    = [3.0];


# 4. Create the JOB structure
job = JOB(thermodynamic_database, solid_solution,
          w_names                                   = w_names,
          w_initial_values                          = w_initial_values,
          w_lower_bounds                            = w_lower_bounds,
          w_upper_bounds                            = w_upper_bounds,
          g0_corr_endmembers                        = g0_corr_endmembers,
          g0_corr_initial_values                    = g0_corr_initial_values,
          g0_corr_lower_bounds                      = g0_corr_lower_bounds,
          g0_corr_upper_bounds                      = g0_corr_upper_bounds,
          algorithm                                 = "NelderMead",          # use "NelderMead" (recommended) or "ParticleSwarm"
          number_iterations_max                     = 5000,                     # maximum number of iterations
          normalization                             = true,                     # normalize the variables to optimise
          number_constraints_max                    = nothing,                  # maximum number of constraints to use for the inversion
          max_time_seconds                          = 180);                     # maximum time in seconds for the inversion

# 6. Check JOB
print_job(job)

# 7. Import the constraints
constraints = read_constraints_from_yml("data/gen_data_FPWMP_biotites.yml");

# filter the constraints to only include the ones, where the phase is biotite is stable
constraints = filter(c -> in("bi",c.assemblage), constraints);

# Call the inversion subroutine
res, norm = inversion(job, constraints, loss_f=chi_squared)

println("----------------------------------------")
println("\nOptimised parameters:\n")
if job.normalization == true
    println(res.minimizer .* norm)
else
    println(res.minimizer)
end
