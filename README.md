
# Simulation study: Augmented-data vs. latent projection

## Description

This is [R](https://www.R-project.org/) code for a simulation study comparing the augmented-data and the latent projection implemented in [**projpred**](https://mc-stan.org/projpred/).

The main file is `sim.R` which is written in a way (e.g., using explicit `print()` statements) so that it can be `source()`-d.
In principle, this file runs the simulation itself and also performs the post-processing.
However, there are some caveats which are described hereafter.
The user options mentioned there can be found in section "User options" of `sim.R`.

For the simulation itself, an initial `brmsfit` is needed, which is created by running once with user option `only_init_fit` set to `TRUE`.
Afterwards, the simulation itself can be run by setting `only_init_fit <- FALSE` and `use_existing_res <- FALSE`.

The post-processing part also includes the creation of plots.
In any case, plots are written to `.pdf` files, but there is also the user option `with_tikz` controlling whether plots should also be written to `.tex` files (by using package [**tikzDevice**](https://CRAN.R-project.org/package=tikzDevice)).
For systems on which no LaTeX distribution is available, `with_tikz` needs to be set to `FALSE`.
By setting `use_existing_res <- TRUE` (and `with_tikz <- TRUE`), the post-processing part can be run again for the same simulation results on a different system featuring a LaTeX distribution.

In summary, the following steps need to be followed if the system on which the simulation itself is run does not have a LaTeX distribution available:

1.  Run `sim.R` with `only_init_fit <- TRUE`.
    User options `use_existing_res` and `with_tikz` are irrelevant in this step.
1.  Run `sim.R` with `only_init_fit <- FALSE`, `use_existing_res <- FALSE`, and `with_tikz <- FALSE`.
1.  On a different system featuring a LaTeX distribution:
    Run `sim.R` with `only_init_fit <- FALSE`, `use_existing_res <- TRUE`, and `with_tikz <- TRUE`.

## Legal information

        simauglat: Simulation study for comparing augmented-data and latent projection in projpred
        Copyright (C) 2022  Frank Weber

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <https://www.gnu.org/licenses/>.
