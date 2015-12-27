# PANGEA.HIV.sim
HIV Phylogenetic Simulation Model for the PANGEA-HIV methods comparison exercise.

### Overview
The Phylogenetics And Networks for Generalised HIV Epidemics in Africa (PANGEA-HIV) consortium 
aims to provide viral sequence data from across sub-Saharan Africa, and to evaluate their viral 
phylogenetic relationship as a marker of recent HIV transmission events [9]. 

Research groups were invited to participate in a blinded methods comparison exercise on simulated 
HIV sequence data sets to test the performance of current phylogenetic methods before their application 
on real data. The methods comparison exercise was announced in October 2014, and closed in August 2015. 
Research groups were asked to estimate 
* reductions in HIV incidence during a typical community-based intervention in sub-Saharan Africa,
* the proportion of HIV transmissions arising from individuals during their first three months of infection (early transmission) at the start of the intervention from sequence data.

Secondary objective were to evaluate the merits of full genome sequence data, the 
impact of various aspects of sequence sampling, and the impact of the proportion of transmissions 
originating from outside the study area.

Generalised HIV-1 epidemics were simulated for a relatively small village population of ~8,000 
individuals and a larger regional population of ~80,000 individuals from two structurally different, agent-based 
epidemiological models. The regional simulation captures individual-level HIV transmission dynamics in a larger 
regional population that is broadly similar to a site (cluster) of the HPTN071/PopART HIV prevention trial in South 
Africa [Hayes]. 

This code repository contains the simulation code to generate sequence data and phylogenetic trees 
that correspond to epidemic scenarios under the regional model. 

### Authors of PANGEA.HIV.sim

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

`PANGEA.HIV.sim` contains just two functions, `sim.regional.args` and `sim.regional`. We now describe these briefly, and the various options in more detail below.
* `sim.regional.args`: Set all input arguments to the simulation. Returns a `data.table` with columns `stat` (name of input argument) and `v` (value of input argument).
* `sim.regional`: Create a UNIX shell script that can be called to generate a single simulation. Takes as input a directory name in which the simulations are to be created, and a `data.table` with all input arguments. Returns the file name to the shell script.

Fire up `R`, and type 

```r
library(PANGEA.HIV.sim)
?sim.regional
```

This will show a first example on how to use the simulation code, followed by code to simulate the PANGEA-HIV 
sequences data sets to address the primary objectives of the methods comparison exercise, and code to simulate 
the PANGEA-HIV tree data sets to address the secondary objectives of the methods comparison exercise.


# Output of the simulation

One simulation produces the following files:

*File name*         | *Description*
--------------------|--------------------
.fa                 | Fasta file of aligned, simulated sequences. One for each gene.
_metadata.csv       | Information on sampled individuals.
_SURVEY.csv         | Cross-sectional surveys conducted on a random subset of the simulated population.
_DATEDTREE.newick   | Trees in newick format. Each tree corresponds to the simulated viral phylogeny among sampled individuals of one simulated transmission chain. One tree per line.

### Available variables in the simulated cross-sectional surveys

*Group variable*    | *Description*
--------------------|--------------------
YR                  | Time at which the survey was taken 
Gender              | Gender
AGE                 | Age

*Count by group*    | *Description*
--------------------|--------------------
ALIVE_N             | Number of sexually active individuals
ALIVE_AND_DIAG_N    | Number of diagnosed individuals
ALIVE_AND_ART_N     | Number of individuals that started ART 
ALIVE_AND_SEQ_N     | Number of individuals with a sequence

### Available variables on simulated individuals
￼
*Individual level variable* | *Description*
-------------------- | ------------------------------------------
IDPOP | Identifier of individual
GENDER | Gender NA if archival sequence

x

# Evolutionary model component

Available full genome subtype C sequences from sub Saharan Africa in the Los Alamos Sequence database were used to seed the simulation of HIV-1C viral variants in the study population. Separate viral introductions occurred at the introduction of HIV into the population in 1980, and thereafter at a constant rate over time. Viral introductions are simulated from a seed sequence. The date of the seed sequence was chosen so that the TMRCA of the simulated sequences falls is consistent with current phylogenetic estimates (Walker, Pybus, Rambaut, & Holmes, 2005).

Within-host lineages that are not part of transmission lineages have a higher evolutionary rate in the simulation (Alizon & Fraser, 2013; Vrancken et al., 2014).

A few ‘archival’ sequences are sampled for the period 1985 to 1999. Since 2000, sequences are more systematically sampled from infected individuals as part of HIV surveillance. Since 2015 until the end of the simulation, the population is more intensely sampled, with roughly the same number of sequences per year. Except of a few simulations, the sequence coverage of the HIV epidemic at the time point 2020.0 is below 10%. Sequence coverage may vary between some data sets. The number of sampled individuals is fixed to 1600 (or 3200 in a few data sets).

For some data sets, the simulated viral phylogenies are provided for each transmission chain with at least one sampled individual in file “*_DATEDTREE.newick”. Each transmission chain phylogeny is in newick format, one phylogeny per line. Tip names match to individual level identifiers in file “*_metadata.csv”. For the remaining data sets, simulated sequences are provided. Sequence names match to individual level identifiers in file “*_metadata.csv”.



