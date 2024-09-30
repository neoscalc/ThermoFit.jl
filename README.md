# ThermoFit

[![Tests](https://github.com/neoscalc/ThermoFit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/neoscalc/ThermoFit.jl/actions/workflows/CI.yml)

ThermoFit is a Julia project for inversion to obtain thermodynamic properties of petrological phases.

This proof-of-concept version of Thermofit.jl was developed for the submission of a research project "An integrated framework for developing accurate and fast observation-driven thermodynamic models" to the Swiss National Science Foundation in October 2024. 

### Installation

```julia-repl
julia>]
  pkg> add https://github.com/neoscalc/ThermoFit.jl
```

### Quick start guide

See in [example](https://github.com/neoscalc/ThermoFit.jl/blob/main/example/minim_w_small.jl) for further information.

```julia
using ThermoFit

CST = global_constants()
PARAMS = global_parameters()

JOB = job(thermodynamic_database,
          solid_solution,
          w_names, w_initial_values, w_lower_bounds, w_upper_bounds,
          algorithm,
          number_iterations_max,
          normalization,
          number_constraints_max,
          max_time_seconds)

inversion_run(JOB, constraints)
```


### Funding

The development of this package was partly supported by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 850530).

ThermoFit was created by Philip Hartmeier & Pierre Lanari; February 2024 (Cassis, France)
