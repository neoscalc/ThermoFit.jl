# ---------------------------------------------------------------------------------------------
#                    ------------------------------------------------
#                   |         THERMOFIT 0.0.1 - EXAMPLE FILE         |
#                    ------------------------------------------------
#
# This file shows an example of how to use the ThermoFit package to adjust Margules parameters.
# Description will be added in the future. For now, please refer to the documentation.
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

# 3. Initialise constants and parameters
CST = global_constants();
PARAMS = global_parameters();

# 4. Initialise the data for the inversion
thermodynamic_database = "mp";
solid_solution = "bi";

# start with the Margules for biotite from White et al. (2014)
# allow only the Margules W(phl,annm) to change within a range of 0-60 kJ/mol
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

w_initial_values = [10  0  0 ;
                    4  0  0 ;
                    10  0  0 ;
                    30  0  0 ;
                    8  0  0 ;
                    9  0  0 ;
                    8  0  0 ;
                    15  0  0 ;
                    32  0  0 ;
                    13.6  0  0 ;
                    6.3  0  0 ;
                    7  0  0 ;
                    24  0  0 ;
                    5.6  0  0 ;
                    8.1  0  0 ;
                    40  0  0 ;
                    1  0  0 ;
                    13  0  0 ;
                    40  0  0 ;
                    30  0  0 ;
                    11.6  0  0];

w_lower_bounds = copy(w_initial_values)
w_lower_bounds[1,1] = -100;
w_upper_bounds = copy(w_initial_values)
w_upper_bounds[1,1] = 100;

# Inversion settings
algorithm = "NelderMead";                # use "NelderMead" (recommended) or "ParticleSwarm"
number_iterations_max = 1000;            # maximum number of iterations
normalization = true;                    # normalize the variables to optimise (Margules)
number_constraints_max = 40;             # maximum number of constraints to use for the inversion
max_time_seconds = 60;                   # maximum time in seconds for the inversion

# 5. Create the JOB structure
JOB = job(thermodynamic_database, solid_solution, w_names, w_initial_values, w_lower_bounds, w_upper_bounds, algorithm, number_iterations_max, normalization, number_constraints_max, max_time_seconds);

# 6. Check JOB
job_check_consistency(JOB)

# 7. Import the constraints
constraints = read_constraints_from_yml("data/gen_data_FPWMP_biotites.yml")

# filter the constraints to only include the ones, where the phase is biotite is stable
constraints = filter(c -> in("bi",c.assemblage), constraints)

# Call the inversion subroutine
res, norm = inversion_run(JOB, constraints)

# println(JOB.w_initial_values) # Note that we update JOB.w_initial_values in this version
job_check_consistency(JOB)
println("Margules parameters:")
if JOB.normalization == true
    println(res.minimizer .* norm)
else
    println(res.minimizer)
end
