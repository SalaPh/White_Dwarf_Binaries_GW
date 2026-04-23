# White_Dwarf_Binaries_GW
Analysis of the gravitational waves (GW) signal from white dwarf binaries (WDBs) detectable with MAGIS Space atom-interferometer. The study focuses on the signal upon merger and can be used to study the future parameter reconstruction of WDB with MAGIS Space.

This code was initially developed as [AIMforGW](https://github.com/sbaum90/AIMforGW.git) by Dr. Sebastian Baum for the project [2309.07952](https://arxiv.org/abs/2309.07952) and then customised for the merging WDB. The reference to this work is [2510.xxxxx](https://arxiv.org/abs/2510.xxxxx).



## Overview of repository:

The repository is divided into the main analysis codes, notebooks, and additional folders containing stored data or supplementary scripts. The notebooks are designed so that they can be used and run to obtain results, without dealing with the other components. 

### Main codes:
- `ParamEstimator_SpaceAI.py` and `ParamEstimator_GroundAI.py`: include the main structure for simulating the GW signature of an evolving binary in frequency space (considering GW evolution only), and can output the Signal-to-Noise Ratio (SNR) and Fisher matrix of the detected signal by a space-based and ground-based interferometer, respectively. It is set to compute this for MAGIS detectors.
- `Parameters_Precision.py` and `Signal_Disappearance.py`: make use of the previous script to compute, respectively, the detectable parameter reconstruction and the alert time for signal disappearance of a series of WDB, which can be arranged in a $\mathcal{M}_c$ - $d_L$ (chirp mass - luminosity distance) or $M_1$ - $M_2$ (the two mass components) grid.

### Main notebooks:
- `example.ipynb`: runs an example run of the code `ParamEstimator_SpaceAI.py` in a structured and easy-to-use way. The input can be simply modified. The results are then saved in the folder `ExampleOutput`. 
- `Plots.ipynb`: runs both  `Parameters_Precision.py` and `Signal_Disappearance.py` and plots some results. More specifically, an $M_1$ - $M_2$ grid and GW parameters can be chosen in the notebook and run for both codes, which store the results in the folder `Precision_Output`. Afterwards, the parameter reconstruction precision for mass ratio $q$, chirp mass $\mathcal{M}_c$, luminosity distance $d_L$ and sky localisation $\Omega$ is plotted. In addition, the alert time of a WDB disappearance is plotted.
- `Orbit_after_RLOF.ipynb`: in contrast with `ParamEstimator_SpaceAI.py`, this notebook computes an alternative WDB orbit evolution with the inclusion of mass transfer. Only the orbit is computed, not the GW signal. In addition, some important properties, such as the merger time and the phase shift due to mass transfer, are computed and plotted, again for an $M_1$ - $M_2$ grid of WDB. The results are stored in the folder `MT_Orbit_after_RLOF`.

### Folders:
- `Analysis_Scripts`: contains additional codes and information needed for the computation. `Instrument_Sensitivity` includes the sensitivity data of MAGIS Space and also the ground-based precursor. `constants.py` includes some constants used overall. `plot_format.py` includes some formatting functions to obtain nice plots.
`antenna_Funs_satellites.py` computes the orbit and direction of the satellite around Earth and Sun.
`antenna_Funs_ground.py` computes the orbit and direction of a ground-based detector with Earth movement.
`waveform_LO.py` and `waveform_PN.py` include some functions of the orbit evolution of a binary, respectively at leading and higher PN orders.
`helper_funs.py` has some functions useful for the main wrappers. Lastly, `WD.py` includes some functions regarding properties and evolution specific to white dwarfs (WDs) and WDBs.
- `ExampleOutput`: represents the storage folder for `example.ipynb`. Currently, it contains the result of a specific run.
- `Precision_Output`: represents the storage folder for `Parameters_Precision.py` and `Signal_Disappearance.py`, or for `Plots.ipynb` if these codes are run through the notebook. Currently, it contains the result of a specific run, for both parameter reconstruction and signal disappearance.
- `MT_Orbit_after_RLOF`: represents the storage folder for `Orbit_after_RLOF.ipynb`. Currently, it contains the result of a specific run.


## Running the code

The code can be installed simply by cloning the repository. 
The additional installation of the package `pycbc`
is required for running the code.

For a first approach, one could run the notebook `example.py` to understand the basic functioning of the main code `ParamEstimator_SpaceAI.py`
Otherwise, the two other notebooks can be easily run as they are to show the results that have already been computed and stored. In addition, one could start modifying the input parameters for new results.


## Citation

If you use this code, please cite the following paper:

```
@article{Sala:2025uqh,
    author = "Sala, Giona and Brandenstein, Chiara and Baum, Sebastian and Graham, Peter W.",
    title = "{Detecting White Dwarf Binary Mergers with Gravitational Waves}",
    eprint = "2510.19913",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "TTK-25-30",
    month = "10",
    year = "2025"
}
```
