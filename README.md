# Fit to 𝐷0 → 𝐾−𝜋+𝑒+𝑒−
## Table of Contents
- [Introduction](#introduction)
- [Theory and experimental aspects](#theory-and-experimental-aspects)
  - [Physics motivation](#physics-motivation)
  - [Challenge](#challenge)
  - [How is the decay simulated and reconstructed?](#how-is-the-decay-simulated-and-reconstructed)
  - [Signal and background sources](#signal-and-background-sources)
- [File structure and implementation](#file-structure-and-implementation)
  - [How is everything linked](#how-is-everything-linked)
  - [Structure](#structure)
- [Example](#example)
  - [Fit results](#fit-results)
- [Requirements](#requirements)
- [Installation (using Docker) and usage](#installation-using-docker-and-usage)
  - [Getting the Software](#getting-the-software)
  - [Building and running the container (with Docker)](#building-and-running-the-container-with-docker)
  - [Output](#output)
## Introduction
This project is a data-analysis project.  
The final aim of is to perform a fit to data of the decay of interest: 𝐷0 → 𝐾−𝜋+𝑒+𝑒−, in order to evaluate the number of signal and background events expected.  
To do so it is necessary to constrain some parameters in the final fit. These parameters are extracted from preliminary fit to simulated samples.  
The simulation of samples for the preliminary study of distirbutions is performed via the external software "RapidSim", which is a fast Monte Carlo simulator, available at the following link:  

```bash
https://github.com/gcowan/RapidSim.git
``` 
However, for our pourpose, RapidSim requires a specific setup (sampling distributions, smearing options, information on the com energy and acceptance, ...) to properly match the configuration and properly reproduce the decay reconstrued at the LHCb detector. 

## Theory and experimental aspects: 
This project focuses on the analysis of the rare decay 𝐷0 → 𝐾−𝜋+𝑒+𝑒−.  
### Physics motivation
This decay is relevant in the ocntext of the study of rare charm decays, specifically 𝐷0 → 𝐾+𝐾−𝑒+𝑒− and 𝐷0 → 𝜋+𝜋−𝑒+𝑒−, since the 𝐷0 → 𝐾−𝜋+𝑒+𝑒− is the normalization channel used as a comparison. From a theoretical perspective, rare semileptonic 4-body charm decays can proceed via both long-distance (LD) and short-distance (SD) interactions. But SD interactions, that proceed via FCNC, are suppressed in the SM by the GIM mechanism, thus they are sensitive to NP physics, making the study of these decay really promising.   
### Challenge
A key experimental challenge arises from the presence of electrons in the final state, which are subject to bremsstrahlung when traversing the detector material. This affects the mass resolution and requires dedicated strategies for photon recovery. To account for this, the analysis is split into Brem 0 (no recovered photons) and Brem 1 (at least one recovered photon). These categories exhibit different background compositions: mis-ID and combinatorial for Brem 0, and partially reconstructed and combinatorial for Brem 1. 
### How is the decay simulated and reconstructed?
The decay 𝐷0 → 𝐾−𝜋+𝑒+𝑒− is reconstructed requiring the 𝐷0 to come from the mother decay D*+ → 𝐷0𝜋+. From the tagging of the soft pion (also called, in this simulatoin, spip = soft pion) it is possible to tag the flavour of the D0 meson.  
The yields are extracted from the maximum likelihood fit to the D0 mass candidate distribution reconstructed as the invariant mass of the 4 daughters particles. 
### Signal and background sources
*Brem 1 category*  
- signal: 𝐷0 → 𝐾−𝜋+𝑒+𝑒−  
- background partially reconstructed: 𝐷0 → 𝐾−𝜋+𝜋-e+nu_e  (1 pion is mis-identified as electron and the neutrino is missed)  

*Brem 0 category*    
- signal: 𝐷0 → 𝐾−𝜋+𝑒+𝑒−  
- background mis-id: 𝐷0 → 𝐾−𝜋+𝜋-𝜋+ (2 pions are mis-identified as electrons)


## File structure and implementation
This repository includes 4 sections, i.e. 4 directories, that correspond to the 4 steps performed for the final simulation:      
1. **Setting**  
    +  `settingMother.cpp`: is responsible for generating the kinematic distributions of the D*+ (Dst) mother particle, which serve as input for the RapidSim simulation. Specifically, it creates 1D histograms of the transverse momentum (pT), pseudorapidity (eta), and azimuthal angle (phi) of the Dst. The phi distribution is taken directly from the Monte Carlo input file `inputFile.root`, while the pT and eta distributions are sampled from analytical PDFs—**Johnson SU** for pT and **Gaussian** for eta—using parameters provided in the external text files `params_Dst_brem0.txt` and `params_Dst_brem1.txt`, corresponding to the two electron bremsstrahlung categories (Brem 0 and Brem 1). The code filters events by Brem category, samples the corresponding distributions, and produces output ROOT files containing the resulting histograms, one for each category. These output files are then used by RapidSim to realistically simulate the Dst production kinematics.

    + `settingDaughters.cpp`: is responsible for generating the smearing inputs required by RapidSim to account for detector resolution effects on the daughter particles of the D*+ (Dst) decay. It processes a ROOT file (`inputFile.root`) containing both true and reconstructed momentum values from official Monte Carlo simulations and calculates the resolution as (P - P_TRUE) / P_TRUE for each daughter. The code treats hadrons and electrons separately, and distinguishes between the Brem 0 and Brem 1 categories based on whether the final-state electrons have undergone bremsstrahlung. For hadrons (spip, Km, pip), the resolution is evaluated in five momentum bins with equal statistics, obtained using quantiles. In each bin, a Gaussian fit is performed to extract the sigma value, and these are compiled into a TGraphErrors object that maps the momentum dependence of the resolution. For electrons (lp, lm), the resolution distributions are stored directly as histograms. The output consists of ROOT files containing these graphs and histograms, one set per Brem category. These files are then used by RapidSim, which applies the smearing to the generated daughter particle momenta by adding a shift whose size is sampled based on the provided resolution inputs.
    
    + `SetInput.cpp`: is used to overwrite the default ROOT input files that RapidSim uses for its simulation and sampling stages, replacing them with updated files generated from our custom distributions and smearing configurations. Depending on the selected Brem category (0 or 1), the code accesses the corresponding output files produced by `settingMother.cpp` and `settingDaughters.cpp`, and copies them into the appropriate directories inside the RapidSim setup. Specifically, it copies the hadron smearing file (`Run3_HadronSmearing.root`) and the electron smearing file (`Run3_ElectronSmearing.root`) into the `rootfiles/smear/` directory, and the D*+ (Dst) mother distribution file (`LHCc14.root`) into the `rootfiles/fonll/` directory. This ensures that RapidSim will use our custom kinematic and smearing inputs when running the simulation.
    /RapidConfigFiles : It's a directory which contains configuration files needed for the simulation to provide the correct smearing to the appropriate particle



2. **Simulation** (with RapidSim)  
    + `simulation.sh`: is a Bash script that runs the RapidSim simulation for different D⁰ decay modes, depending on the selected bremsstrahlung category. When executed with an argument (`0` or `1`), it performs the corresponding simulation setup. In mode `1`, it runs simulations for the `D⁰ → Kπee` and `D⁰ → Kππeν` decays. In mode `0`, it simulates `D⁰ → Kπee` and `D⁰ → Kπππ`. For each decay channel, 50'000 events are generated using RapidSim. The resulting `decay_tree.root` files are stored in their respective directories (`D02Kpiee`, `D02Kpipienu`, or `D02Kpipipi`), ready for further analysis.  

3. **PreliminaryFit**   
    The following files perform the preliminary fit to the reconstructed D0 invariant mass from the simulated data via RapidSim  
    The files are similar but are taking separated since for the Brem 0 category a corrections is necessary.
    + `Brem0/fit.cpp`: performs a fit to the simulated signal and misidentified background distributions for the Brem 0 category. The main goal is to extract and save the parameters of the signal and background probability density functions (PDFs), which will later be used in data fitting. The signal is modeled using a double-sided Crystal Ball function and is fitted to the mass distribution of simulated `D⁰ → Kπee` events produced with RapidSim. The background comes from misidentified `D⁰ → Kπππ` decays, where two pions are mis-identified as electrons. This component is modeled with a Bukin function and includes an efficiency correction based on a 2D histogram (`efficiencies_correct.root`) that accounts for the reconstruction and particle ID effects on the mis-identified pions. Both fits are performed in a defined region of the dilepton invariant mass, and the results are saved in text files along with diagnostic plots showing the fitted distributions and their corresponding pull plots.
 
    + `Brem1/fit.cpp`: performs a fit to the signal and partially reconstructed background components of simulated data in the Brem 1 category. The purpose of this fit is to extract the parameter values of the signal and background probability density functions (PDFs), which will later be used in data fits. For the signal, the code fits a Johnson distribution to the D⁰ mass spectrum from `D⁰ → Kπee` simulated decays. For the background, a Bukin function is fitted to events from partially reconstructed `D⁰ → Kππeν` decays. Both samples are filtered to select a specific dilepton invariant mass window. The results are visualized through mass fits and pull distributions, and the extracted fit parameters are saved to text files along with convergence diagnostics and the covariance matrix.

    INPUT: Rapidsim previous output *decay_tree.root* of the proper directory 
    OUTPUT: *params_signal.txt* and *params_bkg.txt*, i.e. the file containing parameters to provide as input for the final fit to data.

4. **FitToData**  
    The following files perform the final fit to the data.  
    The files are similar but are taking separated since for the 2 categories the fit models are different. 
    + `Brem0/finalFit.cpp`: performs the final fit of the D⁰ mass distribution in the Brem 0 category using simulated events. The total model is an extended sum of three components: the **signal**, modeled with a double-sided Crystal Ball function (`RooCrystalBall`) to capture both core Gaussian resolution and non-Gaussian tails; the **misidentified background**, described by a Bukin function (`RooBukinPdf`) which accounts for asymmetric, peaked background shapes arising from pion misidentification; and the **combinatorial background**, modeled with a first-order Chebychev polynomial (`RooChebychev`). The shape parameters for the Crystal Ball and Bukin components are loaded from previous signal and background fits (`paramsSignal.txt` and `paramsBkg.txt`) and, beside the mean and sigma of the Johnson pdf and the mean peak of the Bukin, they are held constant during the final fit, while the yields and polynomial coefficient are allowed to float. The fit is performed on a dataset stored in `data.root`, with mass values in the range 1600–2150 MeV. The result includes a pull distribution and component-wise visualization, and the full set of fitted yields, correlation matrices, and fit diagnostics are saved for further analysis.

    + `Brem1/finalFit.cpp`: performs the final mass fit for the D⁰ → Kπee simulated sample in the Brem 1 category. The model consists of three components: the **signal** is described by a Johnson distribution (`RooJohnson`), which is well-suited to accommodate asymmetric, non-Gaussian tails often observed in electron final states; the **partially reconstructed background** is modeled using a Bukin function (`RooBukinPdf`), capturing skewed peak-like shapes from decays with missing particles (e.g., neutrinos); and the **combinatorial background** is represented by a first-order Chebychev polynomial (`RooChebychev`). Parameters for the signal and background shapes are imported from previous simulation-based fits (`paramsSignal.txt`, `paramsBkg.txt`) and fixed during this final step (all besides the mean of the signal pdf), while only the component yields and the polynomial coefficient are left free to float. The fit is performed on a `DecayTree_BDTcut` dataset stored in `data.root`, covering the D⁰ mass range from 1600 to 2150 MeV. After fitting, the script evaluates the number of signal and background events within a dynamically defined signal region (mean ± 3σ), calculates their statistical uncertainties and correlations, and produces diagnostic plots including fit overlays and pull distributions. All fit results, including yields, covariance matrices, and convergence status, are saved to text files for further interpretation and documentation.


    INPUT: fit parameters to fix from *params_signal.txt* and *params_bkg.txt* of the correct folder and *data.root*   
    OUTPUT: fit plots, fit results and yields for signal and background. 


**WARNING**: settingDaughters.cpp and settingMotherDst.cpp are executed just one time when the container is built to provide the graphs. It is not necessary to run them several time.   

### How is everything linked:
```bash
./run_pipeline.sh X
``` 
`run_pipeline.sh` is the main script that automates the entire analysis pipeline. It takes a single argument `X` (either `0` or `1`) corresponding to the desired Brem category. The script performs four sequential steps: (1) it sets up the correct input files for RapidSim by running `SetInput.cpp`, ensuring the smearing and kinematic distributions match the selected Brem configuration; (2) it runs the simulation by calling `simulation.sh`, generating 50'000 events per decay channel; (3) it performs the preliminary fits to signal and background components using `fit.cpp` within the appropriate `PreliminaryFit/BremX/` directory; and (4) it executes the final fit using `finalFit.cpp` in `FitToData/BremX/`, combining all modeled components to extract yields and fit diagnostics. This script provides a convenient way to reproduce the full workflow for either Brem 0 or Brem 1 with a single command:

### Structure  
The previously explained structure can be observe here:
```text
.
├── Dockerfile
├── FitToData
│   ├── Brem0
│   │   ├── data.root
│   │   └── finalFit.cpp <---
│   └── Brem1
│       ├── data.root
│       └── finalFit.cpp <---
├── PreliminaryFit
│   ├── Brem0
│   │   ├── efficiencies_correct.root
│   │   └── fit.cpp <---
│   └── Brem1
│       └── fit.cpp <---
├── README.md
├── Setting
│   ├── RapidSimConfigFiles
│   │   ├── Km
│   │   ├── lm
│   │   ├── lp
│   │   ├── pip
│   │   └── spip
│   ├── SetInput.cpp <---
│   ├── inputFile.root
│   ├── params_Dst_brem0.txt
│   ├── params_Dst_brem1.txt
│   ├── settingDaughters.cpp <---
│   └── settingMotherDst.cpp <---
├── Simulation
│   ├── D02Kpiee
│   │   ├── Kpiee.config
│   │   └── Kpiee.decay
│   ├── D02Kpipienu
│   │   ├── Kpipienu.config
│   │   └── Kpipienu.decay
│   ├── D02Kpipipi
│   │   ├── Kpipipi.config
│   │   └── Kpipipi.decay
│   └── simulation.sh <---
└── run_pipeline.sh <---
``` 

***Note***:`<---` = code mentioned

## Requirements  
The project is based on:
- Root framework
- C++
- RooFit
- RapidSim

All encapsulated in:
- Docker (through `Dockerfile`)

## Installation (using Docker) and usage  
### Getting the Software
To get the software `clone` from the official repo and enter in the software folder:

```bash
git clone https://github.com/SofiZala/SoftwareComputingProject.git
cd SoftwareComputingProject
``` 
### Building and running the container (with Docker)
The following steps illustrates how to run properly the simulation and observe the results.  
1. Build the container:  
```bash
docker build -t rapidsim-env .
``` 

2. Run the simualation:  
You can choose which brem fit you want to observe:  
```bash
docker run --name my_container -it rapidsim-env 1
``` 
for the Brem 1 simualtion, or
```bash
docker run --name my_container -it rapidsim-env 0
``` 
for the Brem 0 simulation.  

**WARNING**: During the simulation you may observe some warnings. Don't worry, they come from the RapidSim original code and they do not affect the result of the simulation.  


3. Save the output to look at the plots and results: 
```bash
mkdir -p ./output
``` 
```bash
docker cp my_container:/home/user/FitToData/BremX/FitData.png ./output
``` 
```bash
docker cp my_container:/home/user/FitToData/BremX/FitData_results.txt ./output
``` 
Remember to substitute Brem**X** with 1 or 0! 

**[Optional]**
4. If you want to look at the input information and at the plots provided to RapidSim, you can copy them from their respective directories:  
```bash
mkdir -p ./input
``` 

```bash
docker cp my_container:/home/user/Setting ./input
``` 
```bash
docker cp my_container:/home/user/Simulation ./input
``` 
```bash
docker cp my_container:/home/user/PreliminaryFit/BremX ./input
``` 

**Afterward**
```bash
docker rm my_container
```

***Note***: you may need `sudo` permissions to build, run and remove the container and also to copy files from inside the container to the outside. 


### Output:   
The output is the fit to data, the respective plot, fit results and yield for signal and background events. 





