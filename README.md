# PANGEA.HIV.sim
HIV Phylogenetic Simulation Model for the PANGEA-HIV methods comparison exercise.

## Overview
The Phylogenetics And Networks for Generalised HIV Epidemics in Africa (PANGEA-HIV) consortium aims to provide viral sequence data from across sub-Saharan Africa, and to evaluate their viral phylogenetic relationship as a marker of recent HIV transmission events [9]. 

Research groups were invited to participate in a blinded methods comparison exercise on simulated HIV sequence data sets to test the performance of current phylogenetic methods before their application on real data. The methods comparison exercise was announced in October 2014, and closed in August 2015. Research groups were asked to estimate 

    * reductions in HIV incidence during a typical community-based intervention in sub-Saharan Africa,
    * the proportion of HIV transmissions arising from individuals during their first three months of infection (early transmission) at the start of the intervention from sequence data;

Secondary objective were to evaluate the merits of full genome sequence data, the impact of various aspects of sequence sampling, and the impact of the proportion of transmissions originating from outside the study area.

Generalised HIV-1 epidemics were simulated for a relatively small village population of ~8,000 individuals and a larger regional population of ~80,000 individuals from two structurally different, agent-based epidemiological models. The regional simulation captures individual-level HIV transmission dynamics in a larger regional population that is broadly similar to a site (cluster) of the HPTN071/PopART HIV prevention trial in South Africa [Hayes]. This code repository contains the simulation code to generate sequence data and phylogenetic trees that correspond to epidemic scenarios under the regional model. 

## Authors of PANGEA.HIV.sim

* Oliver Ratmann <oliver.ratmann@imperial.ac.uk>
* Mike Pickles <m.pickles@imperial.ac.uk>
* Anne Cori <a.cori@imperial.ac.uk>
* Christophe Fraser <c.fraser@imperial.ac.uk>

# Installation

The simulation code ships as an `R` package and requires third party code before installation:

* GNU Scientific Library http://www.gnu.org/software/gsl/

After these dependencies are installed, type in `R`:

```r
library(devtools)
install_github("olli0601/PANGEA.HIV.sim")
```

# Using the package

`PANGEA.HIV.sim` contains just two functions, `simulate.regional.args` and `simulate.regional`. We now describe these briefly, and the various options in more detail below.
* `simulate.regional.args`: Set all input arguments to the simulation. Returns a `data.table` with columns `stat` (name of input argument) and `v` (value of input argument).
* `simulate.regional`: Create a UNIX shell script that can be called to generate a single simulation. Takes as input a directory name in which the simulations are to be created, and a `data.table` with all input arguments. Returns the file name to the shell script.

Fire up `R`, and type 

```r
library(PANGEA.HIV.sim)
?simulate.regional
```
