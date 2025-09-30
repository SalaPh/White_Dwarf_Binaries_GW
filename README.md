# White_Dwarf_Binaries_GW
Analysis of the gravitational waves (GW) signal from white dwarf binaries (WDBs) detectable with MAGIS Space atom-interferometer. The study focuses on the signal upon merger and can be used to study the future parameter reconstruction of WDB with MAGIS Space.

This code was initially developed as [AIMforGW](https://github.com/sbaum90/AIMforGW.git) by Dr. Sebastian Baum for the project [2309.07952](https://arxiv.org/abs/2309.07952) and then customised for the merging WDB. The reference to this work is [2510.xxxxx](https://arxiv.org/abs/2510.xxxxx).



## Overview of repository:

The repository is divided into the main analysis codes, notebooks, and additional folders containing stored data or supplementary scripts. The notebooks are designed so that they can be used and run to obtain results, without dealing with the other components. 
In order to run the code, 

### Main codes:
- `ParamEstimator_SpaceAI.py`: includes the main structure for simulating the GW signature of an evolving binary in frequency space (considering GW evolution only), and can output the Signal-to-Noise Ratio (SNR) and Fisher matrix of the detected signal, hence compute the parameter reconstruction.
- `Parameters_Precision.py` and `Signal_Disappearance.py`: make use of the previous script to compute, respectively, the detectable parameter reconstruction and the alert time for signal disappearance of a series of WDB, which can be arranged in a $\mathcal{M}_c$ - $d_L$ (chirp mass - luminosity distance) or $M_1$ - $M_2$ (the two mass components) grid.

### Main notebooks:
- `example.ipynb`: runs an example run of the code `ParamEstimator_SpaceAI.py` in a structured and easy-to-use way. The input can be simply modified. The results are then saved in the folder `ExampleOutput`. 
- `Plots.ipynb`: runs both  `Parameters_Precision.py` and `Signal_Disappearance.py` and plots some results. More specifically, an $M_1$ - $M_2$ grid and GW parameters can be chosen in the notebook and run for both codes, which store the results in the folder `Precision_Output`. Afterwards, the parameter reconstruction precision for mass ratio $q$, chirp mass $\mathcal{M}_c$, luminosity distance $d_L$ and sky localisation $\Omega$ is plotted. In addition, the alert time of a WDB disappearance is plotted.
- `Orbit_after_RLOF.ipynb`: in contrast with `ParamEstimator_SpaceAI.py`, this notebook computes an alternative WDB orbit evolution with the inclusion of mass transfer. Only the orbit is computed, not the GW signal. In addition, some important properties, such as the merger time and the phase shift due to mass transfer, are computed and plotted, again for an $M_1$ - $M_2$ grid of WDB. The results are stored in the folder `MT_Orbit_after_RLOF`.

### Folders:
- `Analysis_Scripts`: contains additional codes and information needed for the computation. Ìnstrument_Sensitivity` includes the sensitivity data of MAGIS Space. `constants.py` includes some constants used overall. `plot_format.py` includes some formatting functions to obtain nice plots.
`antenna_Funs_satellites.py` computes the orbit and direction of the satellite around Earth and Sun.
`waveform_LO.py` and `waveform_PN.py` include some functions of the orbit evolution of a binary, respectively at leading and higher PN orders.
`helper_funs.py` has some functions useful for the main wrappers. Lastly, `WD.py` includes some functions regarding properties and evolution specific to white dwarfs (WDs) and WDBs.
- `ExampleOutput`: represents the storage folder for `example.ipynb`. Currently, it contains the result of a specific run.
- `Precision_Output`: represents the storage folder for `Parameters_Precision.py` and `Signal_Disappearance.py`, or for `Plots.ipynb` if these codes are run through the notebook. Currently, it contains the result of a specific run, for both parameter reconstruction and signal disappearance.
- `MT_Orbit_after_RLOF`: represents the storage folder for `Orbit_after_RLOF.ipynb`. Currently, it contains the result of a specific run.

## Installation

The code can be installed simply by cloning the repository. The repository contains a `WD_GWB.yml` file that can be used to create a `conda` environment that contains all the required packages to run the code in `src` and `post_proc`. This is done by running

```
$ conda env create -f WD_GWB.yml
```

## Running the code

The code can be run from the main directory or `src`.
One should start by activating the `WD_GWB` environment as follows:

```
$ conda activate WD_GWB
```

Afterwards, all three scripts can simply be run by doing

```
$ python Create_z_at_age.py
$ python SeBa_pre_process.py
$ python GWB.py param.ini
```

For `Create_z_at_age.py` one only needs to specify a maximum redshift and the number of interpolation points desired.

For `SeBa_pre_process.py`, one only needs to specify the data paths and whether to save the file. Additional datafolders can be made here, and they should be adapted in the main code.

For `GWB.py`, the settings are specified in a parameter file, where `param.ini` is shown as an example in the repository. As an example, the code takes 7-8 minutes to create the examples `SFH1_50_20_*_example.txt`.

## Clarification on some of the formulas

References to equations are as given in [MT1](/references/master_thesis_Seppe_Staelens.pdf), unless otherwise mentioned. The clarifications here are to correct numerical factors in [MT1](/references/master_thesis_Seppe_Staelens.pdf), that have been adapted for [Staelens, Nelemans 2024](https://www.aanda.org/articles/aa/full_html/2024/03/aa48429-23/aa48429-23.html).

### Expressions for Omega

For the bulk, $\Omega$ is calculated as in (2.45), with a different prefactor 8.10E-9 / $S$ where $S$ is a normalization factor introduced through the SeBA code: SeBa takes as input a certain amount of mass $S$ available for star formation. This factor $S$ is 1.5E6 in [MT1](/references/master_thesis_Seppe_Staelens.pdf), 4E6 in [Staelens, Nelemans 2024](https://www.aanda.org/articles/aa/full_html/2024/03/aa48429-23/aa48429-23.html) and 3.4E6 for [2407.10642](https://arxiv.org/abs/2407.10642).

Similarly, the constant prefactor in (2.48) is changed to 1.28E-8 / $S$, instead of 8.5E-15, for the birth and merger contributions.

### Expressions for z contribution

In the bulk part of the code, the $z$ contributions are saved as $\Omega$ / (some factor), in order not to work with small numbers. This is the reason the $z$ contributions in the merger and birth part are saved as $\Omega$ / (this normalization) as well, to keep the relative contributions the same (as that is all we are interested in). This is admittedly still a bit messy in the code, and should be reworked.

### Expressions for number of binaries

This part is not in [MT1](/references/master_thesis_Seppe_Staelens.pdf). I determined the number of binaries in each $z$-$f$ bin as follows:

$$ N(z, f) = (4 \pi \cdot \chi(z)^2 \cdot \Delta \chi(z)) \cdot n (z, f) , $$

where $n(z, f)$ is the number density of systems in the bin. The latter is given by

$$ n(z, f) = \sum\*k \frac{\psi(z; k)}{S} \cdot \tau(z, f; k) .$$

In this expression, $\psi$ is again the SFH, determined at the birth time of the system and normalized by $S$ as explained above; and $\tau$ is the time it takes the system to traverse the bin. The reasoning is that all systems produced in the past, during a time corresponding to $\tau$, will have moved to the bin under consideration.

In the code, $\tau$ is calculated in Myr, and therefore multiplied by 1E6 as $\psi$ has units of 1/yr.

## Citation

If you use this code, please cite the following paper:

```
@article{
}
```
