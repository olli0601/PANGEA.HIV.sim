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

Secondary objectives were to evaluate the merits of full genome sequence data, the 
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

****

# Output of the simulation

One simulation produces the following files:

*File name*             | *Description*
------------------------|--------------------
.fa                     | Fasta file of aligned, simulated sequences. One for each gene.
_metadata.csv           | Information on sampled individuals.
_SURVEY.csv             | Cross-sectional surveys conducted on a random subset of the simulated population.
_DATEDTREE.newick       | Trees in newick format. Each tree corresponds to the simulated viral phylogeny among sampled individuals of one simulated transmission chain. One tree per line.
_SIMULATED_INTERNAL.R  | Further internal R objects from which the output files are built. Part of the _Internal.zip file.

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

*Individual level variable* | *Description*
-------------------- | ------------------------------------------
IDPOP                | Identifier of individual
GENDER               | Gender (NA if archival sequence)
DOB                 | Date of birth (NA if archival sequence)
DOD                 | Date of death (NA if alive at end of simulation)
DIAG_T              | Time of diagnosis (NA if archival sequence)
DIAG_CD4            | CD4 count at diagnosis (NA if archival sequence)
ART1_T              | ART start date (NA if ART not started)
ART1_CD4            | CD4 count at ART start (NA if ART not started)
TIME_SEQ            | Date sequence taken
RECENT_TR           | Y if transmission occurred at most 6 months after diagnosis N otherwise

****

# Evolutionary model component

The model simulates viral phylogenies of sampled transmission chains in the regional population, as well as
HIV sequences comprising the gag gene (from p17 start, length 1440 nt), the pol gene (from protease start, length 2844) 
and the env gene (from TVA signal peptide start, length 2523, V loops excluded) of sampled individuals.

The model comprises the following components, which are described in detail below:

1. [Viral introductions and seed sequences](#viral-introductions-and-seed-sequences)

2. [Dated viral phylogenies](#dated-viral-phylogenies)

3. [Evolutionary rate model](#evolutionary-rate-model)

4. [Substitution model](#substitution-model)

5. [Sequence sampling model](#sequence-sampling-model)

## Viral introductions and seed sequences

#### Viral introductions

The evolutionary model starts in 1980. By 1980, between 100-200 infected individuals exist in the epidemic simulation. 
The evolutionary model considers these individuals as index cases. 
After 1980, viral introductions occur in proportion to the number of new infections per year (parameter `epi.import`). 
New index cases representing viral introductions into the regional population are selected at random among existing infected individuals. 

By default, we set `epi.import=0.05` in line with baseline assumptions of the HPTN071/PopART model 
(Cori A, Ayles H, Beyers N, Schaap A, Floyd S, et al. PLoS One 2014). Optionally, we also used `epi.import=0.2` reflecting
estimates of substantial viral introductions from outside the Rakai community cohort in Uganda (Grabowski MK, Lessler J, Redd AD, et. al. PLoS Med 2014).
Sources of index cases are coded with a negative population ID `IDPOP`. The corresponding transmission chains
can be reconstructed from data.table `df.trms` (in file ending `_SIMULATED_INTERNAL.R`). Column `IMPORTp` of data.table `df.epi` 
(in file ending `_SIMULATED_INTERNAL.R`) lists the simulated proportion of viral introductions per year.

#### Seed sequences

Each index case is assigned a seed sequence. The date associated with the sequence is either 1980.0 or, for viral introductions 
after 1980, the time of infection of the index case. To obtain a well-defined tree, the seed sequences are evolved from one 
starting sequence (parameter `startseq.mode`). The starting sequence is randomly sampled from a pool of pre-specified sequences
by date (parameter `index.starttime.mode`, see also next section). To begin constructing the simulated phylogeny, 
the starting sequence is connected to each index case by a root edge. The length of the root edge is the time between 
the date of the starting sequence and the time of viral introduction.

By default, one starting sequence is sampled at random from pre-specified starting sequences from 1970 (`startseq.mode="one"`, 
`index.starttime.mode="fix1970"`). The date of the starting sequence is chosen so that the estimated TMRCA of the simulated sequences 
under a root-to-tip divergence plot is consistent with current phylogenetic estimates of the origin of subtype C sequences in 
sub-Saharan Africa (Walker, Pybus, Rambaut, & Holmes, 2005).

The following figure shows part of a simulated phylogeny under default parameters: a single starting sequence evolves into four sub-clades. 
These four sub-clades correspond to four transmission chains after their introduction in the regional population. 
The four sub-clades are connected with a root edge from the index case to the starting sequence. Under default parameters, it is by design 
easy to correctly identify sub-clades that correspond to distinct transmission chains. 

![alt tag](https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_viralintro.png)

#### Starting sequences

The starting sequence was selected from a pool of pre-specified, historical full genome subtype C sequences. This pool was generated
through ancestral state reconstruction with BEAST 1.8 from 390 full genome subtype C sequences from the Los Alamos Sequence database. The alignment
of these 390 sequences was manually curated and is [available here](https://github.com/PangeaHIV/HPTN071sim/tree/master/raw_rootseq). 
The pool of dated historical sequences is [available here](https://github.com/olli0601/PANGEA.HIV.sim/blob/master/inst/misc/PANGEA_SSAfgBwhRc-_140907_n390_AncSeq.R). The following figure shows a histogram of the available starting sequences by their sampling date.

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_PANGEA_SSAfgBwhRc-_140811_n390_geneENV_AncSeq_Times.png" width="66%">

## Dated viral phylogenies

The topology of viral phylogenies does not necessarily correspond to the transmission tree, especially when viral infections persist 
life-long (Pybus OG, Rambaut A, Nat Gen Rev 2009). To allow for such incongruencies, we used a particular within- and between host 
coalescent model that is more fully described elsewhere (Didelot X, Gardy J, Colijn C, Mol Bio Evol 2014; Hall, MD PhD thesis 2015). 

Briefly, continuing with every index case, viral phylogenies with branch lengths in calendar time are generated through recursive 
application of a neutral within-host coalescent model. The infection time of the index case is considered as root of the 
within-host phylogeny of the index case, and any onward transmission events or sampling events as tips. 
Infection and onward transmission times are taken from the epidemic simulation, and sampling times are described further below. 
Under these tip and date constrains, the within-host phylogeny of the index case is simulated assuming an increasing
effective population size (lognormal model, parameters `v.N0tau`, `v.r`, `v.T50`). For each new infection, the process is repeated and the 
within-host phylogenies of newly infected individuals are concatenated to the corresponding transmission tips of their transmitter. 
The model assumes that a single transmitted virion leads to clinical infection of
the newly infected individual. For each index case, the simulation produces a dated viral phylogeny that is rooted at the index case and has as tips
the sampling times of all individuals in the same transmission chain that are sampled.

The following figure shows the lognormal model of the within-host effective population size under default 
parameters `v.N0tau=1`, `v.r=2.851904`, `v.T50=-2`. These were chosen such that the initial population size is 1 and that 
the final effective population size is broadly similar to estimates typically obtained with a BEAST Skyline model.

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_EffPopSize.png" width="66%">

Source code of the VirusTreeSimulator was provided by Matthew Hall. 
The lognormal effective population size model is inheried from BEAST, `BEAST::LogisticGrowthN0`.
Further details are [available here](https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/VirusTreeSimulator_maths.pdf), 
with time expressed as negative time since infection. 


## Evolutionary rate model

To convert branch lengths from calendar times into substitutions per site, evolutionary rates are assigned to each branch. 
Non-transmission lineages are assigned higher evolutionary rates than transmission lineages (Alizon & Fraser, 2013; Vrancken et al., 2014).
Evolutionary rates of non-transmission lineages are sampled from a log normal model with mean `wher.mu` and 
standard deviation `wher.sigma` on the log scale. Evolutionary rates of transmission lineages are sampled from a 
log normal model with mean `bwerm.mu` and standard deviation `bwerm.sigma` on the log scale. The following figure shows the
evolutionary rate model for transmission and non-transmission lineages under default parameters `wher.mu=log(0.00447743)-0.5^2/2`, `wher.sigma=0.5`, 
`bwerm.mu=log(0.002239075)-0.3^2/2`, `bwerm.sigma=0.3`. These were chosen such that XXX

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_ERmodel.png" width="66%">

Optionally,

Each transmission chain phylogeny has a long root branch dating to between 1945-1960 in order to select a starting sequence. The evolutionary rate for these root branches was set to the estimated overall mean evolutionary rate across all genes (2.2e-3 subst/site/year).

Tip branches in each transmission chain phylogeny are modelled to evolve under a high within host rate. All other branches are part of transmission lineages and are modelled to evolve under a slower rate. Evolutionary rates are drawn independently for each individual from lognormal distributions with means 2.2e-3 subst/site/year (along transmission lineages) and 4.4e-3 subst/site/year (along tip branches). 32% of all branches are tip branches. This is specified with wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13 in the rPANGEAHIVsim.pipeline.args function.



## Substitution model

## Sequence sampling model


A few ‘archival’ sequences are sampled for the period 1985 to 1999. Since 2000, sequences are more systematically sampled from infected individuals as part of HIV surveillance. Since 2015 until the end of the simulation, the population is more intensely sampled, with roughly the same number of sequences per year. Except of a few simulations, the sequence coverage of the HIV epidemic at the time point 2020.0 is below 10%. Sequence coverage may vary between some data sets. The number of sampled individuals is fixed to 1600 (or 3200 in a few data sets).

For some data sets, the simulated viral phylogenies are provided for each transmission chain with at least one sampled individual in file “*_DATEDTREE.newick”. Each transmission chain phylogeny is in newick format, one phylogeny per line. Tip names match to individual level identifiers in file “*_metadata.csv”. For the remaining data sets, simulated sequences are provided. Sequence names match to individual level identifiers in file “*_metadata.csv”.



