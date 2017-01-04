# PANGEA.HIV.sim
HIV Phylogenetic Simulation Model for the PANGEA-HIV methods comparison exercise.

#### Overview
The Phylogenetics And Networks for Generalised HIV Epidemics in Africa (PANGEA-HIV) consortium 
aims to provide viral sequence data from across sub-Saharan Africa, and to evaluate their viral 
phylogenetic relationship as a marker of recent HIV transmission events 
(Pillay D, Herbeck J, Cohen MS, et al. Lancet Infect Dis 2015). 

Research groups were invited to participate in a blinded methods comparison exercise on simulated 
HIV sequence data sets to test the performance of current phylogenetic methods before their application 
on real data (Ratmann O, Hodcroft E, Pickles M, et al. Molecular Biology Evolution 2016). The methods comparison exercise was announced in October 2014, and closed in August 2015. 
Research groups were asked to estimate 
* reductions in HIV incidence during a typical community-based intervention in sub-Saharan Africa,
* the proportion of HIV transmissions arising from individuals during their first three months of infection (early transmission) at the start of the intervention from sequence data.

Secondary objectives were to evaluate the merits of full genome sequence data, the 
impact of various aspects of sequence sampling, and the impact of the proportion of transmissions 
originating from outside the study area.

Generalised HIV-1 epidemics were simulated for a relatively small village population of ~8,000 
individuals and a larger regional population of ~80,000 individuals from two structurally different, agent-based 
epidemiological models. The regional simulation captures individual-level HIV transmission dynamics in a 
regional population that is broadly similar to a site (cluster) of the HPTN071 (PopART) HIV prevention trial in South 
Africa (Hayes R, Ayles H, Beyers N, et al., Trials 2014).

This repository implements the regional model to generate sequence data and phylogenetic trees 
for a variety of generalised HIV epidemics and intervention scenarios. It is based on the HPTN071/PopART individual based model, see also here: https://github.com/mrehp2/POPART-IBM/tree/PangeaSim.

#### Authors of PANGEA.HIV.sim

* Oliver Ratmann <oliver.ratmann@imperial.ac.uk> (maintainer)
* Mike Pickles <m.pickles@imperial.ac.uk>
* Anne Cori <a.cori@imperial.ac.uk>
* Christophe Fraser <c.fraser@imperial.ac.uk>
 
****

# Using the package

#### Installation

The simulation code ships as an `R` package and requires third party code before installation:

* GNU Scientific Library http://www.gnu.org/software/gsl/

After these dependencies are installed, type in `R`:

```r
library(devtools)
install_github("olli0601/PANGEA.HIV.sim")
```

#### Functionality

`PANGEA.HIV.sim` contains just two functions, `sim.regional.args` and `sim.regional`. We now describe these briefly, and the various options in more detail below.
* `sim.regional.args`: Set all input arguments to the simulation. Returns a `data.table` with columns `stat` (name of input argument) and `v` (value of input argument).
* `sim.regional`: Create a UNIX shell script that can be called on the command line to generate a single simulation. Takes as input a directory name in which the simulations are to be created, and a `data.table` with all input arguments. Returns the file name to the UNIX shell script. 

The package was developped on MacOS X and should also work on LINUX. Windows is not supported, and support is not planned.

Fire up `R`, and type 

```r
library(PANGEA.HIV.sim)
?sim.regional
```

This will show a first example on how to use the simulation code, followed by code to simulate the PANGEA-HIV 
sequences data sets to address the primary objectives of the methods comparison exercise, and code to simulate 
the PANGEA-HIV tree data sets to address the secondary objectives of the methods comparison exercise. 
To re-cap the first example:

```r
library(PANGEA.HIV.sim)
outdir          <- getwd()                                          #set to a new empty directory; dir name must not contain whitespace, brackets, etc
pipeline.args	<- sim.regional.args( 	seed=42,                    #random number seed for reproducibility
                                        yr.end=2020,				#end of simulation
                                        s.PREV.max.n=1600,          #number of sequences
                                        s.INTERVENTION.prop=0.5,    #proportion of sampled sequences after intervention start in 2015
                                        epi.acute='high',           #frequency of early infections (high or low)
                                        epi.intervention='fast',    #intervention scale-up (none, slow or high)
                                        epi.import=0.05 )			#proportion of transmissions from outside the regional population
cat(sim.regional(outdir, pipeline.args=pipeline.args))              #run this script from the command line
```

#### Output of the simulation

###### Output files

*File name*             | *Description*
------------------------|--------------------
.fa                     | Fasta file of aligned, simulated sequences. One for each gene.
_metadata.csv           | Information on sampled individuals.
_SURVEY.csv             | Cross-sectional surveys conducted on a random subset of the simulated population.
_DATEDTREE.newick       | Trees in newick format. Each tree corresponds to the simulated viral phylogeny among sampled individuals of one simulated transmission chain. One tree per line.
_SIMULATED_INTERNAL.R  | Further internal R objects from which the output files are built. Part of the _Internal.zip file.

###### Available variables in the simulated cross-sectional surveys

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

###### Available variables on simulated individuals

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

#### Simulations used in the PANGEA methods comparison exercise

The 33 simulated data sets from the village and regional models are [available here in folder 201502](http://bit.ly/PANGEAHIVsim).

****

# Epidemiological model component

## Regional model

The regional simulations were produced using the individual-based HPTN071 (PopART) model version 1.1. This is a stochastic individual-based model of a heterosexual population of approximately 80,000 individuals aged 13 years and older. At a given time in the simulation a combination HIV prevention intervention is initiated. This intervention is a combination HIV prevention intervention, qualitatively similar to the HPTN071/PopART trial, although coverage levels used were not informed by trial data. The model is calibrated to have a similar HIV prevalence to that observed in the South African sites of the trial. The outputs from the model include the full HIV transmission chain within the simulated community, which is then passed to the evolutionary model component **(ADD LINK)** to generate phylogenies and viral sequence data sets.

The model comprises the following components, which are described in detail below:

 * Demographics
 * Sexual partnerships and HIV infection
 * Disease progression
 * HIV testing
 * Male circumcision
 * Antiretroviral therapy

### Demographics

Individuals enter the modelled population continuously at age 13, leaving the population when they die, either from HIV- or non HIV-related reasons. The overall population is slowly growing in size. 

On entering the model, each individual is assigned a set of characteristics including gender, a sexual risk group (there are three risk groups, corresponding to different allowed maximum numbers of concurrent partnerships), and (for men only) if they are already circumcised (13% assumed to be circumcised before the age of 13). All individuals are assumed HIV-negative when entering the simulation. 

### Sexual partnerships and HIV infection
Once they have entered the model, they may form and break up sexual partnerships. Individuals form partnerships assortatively by risk group, with 10% of one’s partnerships made within a risk group, and the remainder made at random in any group. Partnership formation is also strongly assortative by age; the age mixing matrix was determined from HPTN071/PopART cohort baseline survey data.
HIV-negative individuals are exposed to risk of infection when they are in a serodiscordant partnership. Infection can occur at any time during a serodiscordant partnership, with the risk of infection depending upon the HIV stage of the infected partner, and whether they are on ART or not at that time. For male HIV negative individuals risk of infection also depends on their circumcision status.  

### Disease progression

A set point viral load (SPVL) category (< 4, 4-4.5, 4.5-5, ≥ 5 log<sub>10</sub> copies/μL) is randomly attributed to each newly infected individual. A newly infected individual firstly undergoes a brief period of early infection (lasting on average 3 months), with higher infectivity. They then enter one of 4 CD4 categories (CD4 >500, 350-500, 200-350, ≤200 cells/mm<sup>3</sup>) according to the distribution from [Cori et al. AIDS 2015][Cori2015] and depending on their SPVL category. Once in that CD4 category they then progress sequentially to lower CD4 categories, at a rate determined by the SPVL category as well as the current CD4 category, before dying from AIDS-related illness as in [Cori et al. AIDS 2015][Cori2015], unless they test for HIV and subsequently initiate ART. 

### HIV testing

HIV testing is divided into two separate rates. Firstly there is a background (standard of care) rate that increases over time, mimicking the historical scale-up of testing in sub-Saharan Africa. Secondly, starting in 2015, there is an intensive annual testing process, representing the HIV-testing component of the combination HIV prevention intervention, where a percentage of the population each year was tested for HIV. 

### Male circumcision

Men testing HIV-negative, who were not previously circumcised, can be circumcised either as part of the background (non-intervention) or intervention process. Once circumcised, susceptibility to HIV is reduced.

### Antiretroviral therapy

Individuals only start ART after a positive HIV test result, although they may be lost to follow-up before this occurs. Once on ART there is a period of early ART (lasting on average 6 months) where the individual remains virally unsuppressed, following which they may either become virally suppressed, or remain on ART but not fully suppressed. Infectivity is reduced when an individual is on ART, but it is more substantially reduced if they are virally suppressed.
Individuals on ART may drop out from treatment. Once outside treatment, individuals may restart treatment once their CD4 goes below 200 cells/mm<sup>3</sup>. 

### Model simulations

The simulation was initially run for 20 years without HIV, to allow partnerships to reach a steady state. In 1980, 0.5% of low-risk and 1% of medium and high-risk individuals are seeded HIV-positive, with HIV transmission occurring from that point onwards. ART is taken to be available from 2004 onwards. The combination prevention intervention runs from 2015 for three years. 20 simulations were run using this model.

### Model calibration

The 20 simulations run using this model differed by uptake and coverage of the combination HIV prevention intervention, or by increased infectivity from early infection. Each run was calibrated to give an incidence close to 2%/year and an HIV prevalence of approximately 15-20% at the start of the intervention. For a given specified set of parameters, calibration was achieved by varying the baseline transmission probability (the probability that an individual with CD4>500 cells/mm<sup>3</sup>, and not on ART, transmits to their partner until the HIV prevalence in 2015 was within the specified range, and HIV prevalence was not declining sharply prior to the intervention.

| Parameter  | Value/distribution used | Reference  |
| ------------- |-------------| -----|
| Timestep	     | 1/48 year	    | Assumption |
| Probability of transmission per timestep for individual with CD4>500cells/mm<sup>3</sup> not on ART (baseline transmission probability)     | Varied to calibrate model incidence/prevalence | Calibrated |
| Relative increase in transmission probability during early transmission phase      | 6.0 (when ~10% of transmissions early); 26.0 (when ~40% of transmissions early)      |   Informed by [Hollingsworth et al. JID 2008][Hollingsworth2008], [Boily et al. LID 2009][Boily2009] and calibrated to get 10% and 40% of transmissions early |
| Relative increase in transmission probability (compared to baseline transmission probability) when CD4 350-500 / 200-350 / ≤200 cells/mm<sup>3</sup> | 1.0 / 1.9 / 3.0      |  Assumption   |
| Effectiveness of ART when virally suppressed	| 0.9	| [Cori et al. PLoS One 2014][Cori2014] |
| Effectiveness of ART when no virally suppressed, or during early ART	| 0.45	| Assumed to be half that when virally suppressed |
| Relative reduction in susceptibility when circumcised	| 0.6	| [Auvert et al. PLoS Med 2005][Auvert2005], [Bailey et al. Lancet 2007][Bailey2007], [Gray et al. Lancet 2007][Gray2007] |
| Duration of early transmission phase (months)	| Uniform(1,5)	| Informed by [Hollingsworth et al. JID 2008][Hollingsworth2008] |
| Duration of CD4 stages  >500, 350-500, 200-350 and ≤200 cells/mm<sup>3</sup>	| Sampled uniformly from ranges from early draft of Cori et al. AIDS 2015 | [Cori et al. AIDS 2015][Cori2015] |
| Relative rate of progression when on ART but virally unsuppressed	| 0.5	| [Cori et al. PLoS One 2014][Cori2014] |
| Relative rate of progression when on ART and virally unsuppressed	| 0 (no progression)	| Assumption |
| Proportion in SPVL group < 4, 4-4.5, 4.5-5, ≥ 5 log<sub>10</sub> copies/μL after seroconversion	| 0.25 in each group	| Informed by [Cori et al. AIDS 2015][Cori2015] |
| Proportion in high / medium / low sexual risk groups when entering population	| Men: 0.5 / 0.4 / 0.1; Women: 0.6 / 0.3 / 0.1	| Assumption |
| Time to partnership dissolution	| Gamma distribution with shape parameter 10 and scale parameter 2.5	| Assumption |

[Hollingsworth2008]: http://dx.doi.org/10.1086/590501
[Cori2014]: http://dx.doi.org/10.1371/journal.pone.0084511
[Cori2015]: http://journals.lww.com/aidsonline/Fulltext/2015/11280/CD4__cell_dynamics_in_untreated_HIV_1_infection__.7.aspx
[Boily2009]: http://dx.doi.org/10.1016/S1473-3099(09)70021-0
[Auvert2005]: http://dx.doi.org/10.1371/journal.pmed.0020298
[Bailey2007]: http://dx.doi.org/10.1016/S0140-6736(07)60312-2
[Gray2007]: http://dx.doi.org/10.1016/S0140-6736(07)60313-4


****

# Evolutionary model component

The model simulates viral phylogenies of sampled transmission chains in the regional population, as well as
HIV sequences comprising the gag gene (from p17 start, length 1440 nt), the pol gene (from protease start, length 2844) 
and the env gene (from TVA signal peptide start, length 2523, V loops excluded) of sampled individuals. For simplicity,
we refer to the concatenated gag, pol, env sequences as full genome sequences.

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

By default, we set `epi.import=0.05` in line with baseline assumptions of the pre-PopART model 
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
sub-Saharan Africa (Walker PR, Pybus OG, Rambaut A, Holmes EC, Infect Gen Evol 2005). These parameters were not changed during the PANGEA simulations.

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

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_StartingSequenceByDate.png" width="66%">

## Dated viral phylogenies

The topology of viral phylogenies does not necessarily correspond to the transmission tree, especially when viral infections persist 
life-long (Pybus OG, Rambaut A, Nat Gen Rev 2009). To allow for such incongruencies, we used a particular within- and between host 
coalescent model that is more fully described elsewhere (Didelot X, Gardy J, Colijn C, Mol Bio Evol 2014; Hall, MD PhD thesis 2015). 

Briefly, continuing with every index case, viral phylogenies with branch lengths in calendar time are generated through recursive 
application of a neutral within-host coalescent model. The infection time of the index case is considered as root of the 
within-host phylogeny of the index case, and any onward transmission events or sampling events as tips. 
Infection and onward transmission times are taken from the epidemic simulation, and sampling times are described further below. 
Under these tip and date constrains, the within-host phylogeny of the index case is simulated assuming an increasing
effective population size (logistic model, parameters `v.N0tau`, `v.r`, `v.T50`). For each new infection, the process is repeated and the 
within-host phylogenies of newly infected individuals are concatenated to the corresponding transmission tips of their transmitter. 
The model assumes that a single transmitted virion leads to clinical infection of
the newly infected individual. For each index case, the simulation produces a dated viral phylogeny that is rooted at the index case and has as tips
the sampling times of all individuals in the same transmission chain that are sampled. 
These phylogenies are all connected via the starting sequence as explained in section 
[Viral introductions and seed sequences](#viral-introductions-and-seed-sequences), and printed to the file ending in `_DATEDTREE.newick`. 

The following figure shows the logistic model of the within-host effective population size under default 
parameters `v.N0tau=1`, `v.r=2.851904`, `v.T50=-2`. These were chosen such that the initial population size is 1 and that 
the final effective population size is broadly similar to estimates typically obtained with a BEAST Skyline model.
These parameters were not changed during the PANGEA simulations.

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_EffPopSize.png" width="66%">

Source code of the VirusTreeSimulator was provided by Matthew Hall. 
The logistic effective population size model is inheried from BEAST, `BEAST::LogisticGrowthN0`.
Further details are [available here](https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/VirusTreeSimulator_maths.pdf), 
with time expressed as negative time since infection. 


## Evolutionary rate model

#### Overall evolutionary rates of transmission and non-transmission lineages

Following the previous steps, the evolutionary model generates a dated viral phylogeny comprising sequences 
of all sampled individuals in the regional population. To evolve sequences from the starting sequence along this phylogeny, 
the branches need to be transformed from calendar times into substitutions per site.
For this purpose, evolutionary rates are assigned to each branch. 
Non-transmission lineages are assigned higher evolutionary rates than transmission lineages 
(Vrancken B, Rambaut A, Suchard MA, et al., PLoS Comput Biol 2014). 

Evolutionary rates of non-transmission lineages are sampled from a log normal model with mean `wher.mu` and 
standard deviation `wher.sigma` on the log scale. Evolutionary rates of transmission lineages are sampled from a 
log normal model with mean `bwerm.mu` and standard deviation `bwerm.sigma` on the log scale. To reduce rate 
variation, it is also -optionally- possible to set the evolutionary rate of the root edges below index cases to the mean
evolutionary rate of transmission lineages with `root.edge.fixed=1`.

By default, we used a mean evolutionary rate of `0.0022` that was estimated as a by-product of the BEAST run on 390 full genome subtype C sequences 
for the pol gene (see section [Starting sequences](#starting-sequences)). 
The mean evolutionary rate of non-transmission lineages was set to twice this value 
(Alizon S, Fraser C, J Retrov 2013; Vrancken B, Rambaut A, Suchard MA, et al., PLoS Comput Biol 2014). The standard deviations of the two log normal 
rate sampling model were chosen so that sequence evolution was relatively clock-like. Specifically, the standard deviations were calibrated so that 
a linear clock model explained about 30-40% of the variation in root-to-tip divergence. The default model parameters are 
`wher.mu=log(0.0044)-0.5^2/2`, `wher.sigma=0.5`, `bwerm.mu=log(0.0022)-0.3^2/2`, `bwerm.sigma=0.3` on the log scale. Note
that the non-transmission lineages are the dead-end sub-clades within an individual that do not lead to onward transmission, but 
along which a sequence was sampled. Since in the model only one sequence is obtained per individual, the non-transmission lineages are a
certain part of the branch just below a tip. To calculate that part, the dated viral phylogeny in section 
[Dated viral phylogenies](#dated-viral-phylogenies) is internally annotated with infection times along each branch.
The following figure illustrates the resulting sampling model of evolutionary rates for transmission and non-transmission lineages under 
the default parameters. These parameters were not changed during the PANGEA simulations.

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_ERmodel.png" width="66%">


#### Relative evolutionary rates by gene and codon position

The above default parameterisation corresponds to the overall rate of evolution of the pol gene.
Relative rate multipliers were used to account for rate differences along codon position and the HIV gag, 
pol and env genes (Alizon S, Fraser C, J Retrov 2013).

The default rate multipliers were obtained as a by-product of the BEAST run on 390 full genome subtype C sequences 
(see section [Starting sequences](#starting-sequences)). The following boxplots show the resulting evolutionary rates along
transmission lineages, stratified by gene and codon position, in a sample simulation.
These parameters were not changed during the PANGEA simulations.

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_geneEvolutionaryRates.png" width="66%">

## Substitution model

A GTR+Gamma substitution model was used to evolve the viral sequences along the viral phylogeny from the starting sequence. 
Each codon and gene were considered separately.
The average number of substitution per branch for each codon and gene were obtained by multiplying the branches in calendar time 
by the respective evolutionary rates in section [Evolutionary rate model](#evolutionary-rate-model).

The rate substitution parameters were obtained as a by-product of the BEAST run on 390 full genome subtype C sequences 
for the pol gene (see section [Starting sequences](#starting-sequences)), and are shown in the figure below. 
[SeqGen/1.3.3](https://github.com/rambaut/Seq-Gen) was used to perform the sequence simulation.
These parameters were not changed during the PANGEA simulations.

<img src="https://github.com/olli0601/PANGEA.HIV.sim/blob/master/man/fig_GTR.png" width="66%">

## Sequence sampling model

#### Sampling of historical sequences

A number of historical sequences can be sampled at random from individuals that are infected during the period 1985 to 1999 
(parameter `s.ARCHIVAL.n`). By default, we set `s.ARCHIVAL.n=50`. This parameter was not changed during the PANGEA simulations.

#### Sampling of sequences since 2000

Since 2000, sequences were randomly sampled at time of diagnosis in proportion to the number of annual new diagnoses.
The proportions of individuals sampled between 2000-2014 and 2015-2020 may differ from each other, and are controlled
by two parameters. The first parameter is the total number of sequences sampled, `s.PREV.max.n`. The second parameter
is the proportion of sampled sequences that are obtained after intervention start in 2015, `s.INTERVENTION.prop`. 
The sampling duration after start of the intervention can be set with `yr.end`.
Column `s.nTOTAL` of data.table `df.sample` (in file ending `_SIMULATED_INTERNAL.R`) lists the number of sampled sequences per year, 
and column `SCOV` in the same data.table lists the sequence coverage among infected and alive individuals by year.
Note this model assumes that all individuals are equally likely sampled.

By default, we set `s.PREV.max.n=1600`, `s.INTERVENTION.prop=0.5` and `yr.end=2020`, which yields a sampling coverage of 8-10% by 
the end of the simulation. With these parameters, sequence coverage increases fast after the start of the intervention.
Optionally, we also used `s.PREV.max.n=1280` (in conjunction with a sampling duration of 3 years), 
`s.PREV.max.n=3200` and `s.INTERVENTION.prop=0.85`. 





