# Where is what to find
This SI contains instruction videos, installation tutorials als examples refinement for the analysis of wide-angle X-ray/neutron (WAXS/WANS) scattering data of non-graphitic carbons (NGCs) using the Ruland & Smarsly (2002) algorithm.

## SI-files
### Supporting information part A – Octave installation and examples
* S1. Overview of all used parameters
* S2. Implementation and calculation time consumption of nu
* S3. Results of the refined samples for the verification of the used software including the calculation times for the WAXS refinements
* S4. Tests for fitting the (004)-region of the LSPP-1200 WANS-data
* S5. Download and usage of iObs
* S6. Installation and updates of Octave
* S7. Usage of Octave
* S8. Example refinement
* S8.9. Common errors during the refinement (see the file Octave/readme.PDF).
* S8.10. Initial refinement script lines 1 – 134
* S8.11. Refinement script lines 1 – 134 after the manual fitting
* S9. Octave cannot always calculate parameter errors - what to do

### Supporting information part B – Correction/fine treatment of WAXS/WANS data and mathematical background
* S11. Incoherent scattering and correction terms
* S12. General intensity correction terms
* S13. Atomic form factors
* S14. Incoherent scattering - theoretical and calculated data
* S15. Background correction for wide-angle neutron scattering

## Example refinement
Hence, the refinement scripts for *Octave* as well as the corresponding "iObs.oct" file can be found. Normally, "WAXS Fit-Routine-IUCr.m" should be used and only if the error bars from the refinement parameters are from interest, but cannot be calculated, try the script "WAXS Fit-Routine-IUCr with additional step all-without-normalization.m".

### Good fit
An example for a "good" *Octave* fit with the important fit regions marked in green (important), yellow (nice to have) and red (unimportant) (Figure 06 of the main article).

### WAXS Steps
An example refinement with given scattering data and examples/scripts for each refinement step.

## Instruction Videos
As the name says, some instruction videos for the installation of *Octave*, fitting of WAXS/WANS data of NGCs and the compilation of an iObs.oct-file.

## *Octave*
The directory *Octave* contains instructions for installing *GNU Octave* on Windows, macOS and Linux. In addition, a refinment script for get the microstructure date of non-graphitic carbons (NGCs) from wide-angle X-ray/neutron (WAXS/WANS) scattering data using the Ruland & Smarsly (2002) algorithm. Also, an instruction for compiling the needed *iObs.oct* file as well as some already compiled *iObs.oct* files are also provided, which can be used to call *iObs* from Octave. For optimal functionality, the file should be renamed to *iObs.oct* after downloading/compiling.
If you never worked with* Octave* so far, it is highly recommended to take a look on an instruction video. Due to the high distribution and availability of *Octave*, there are a lot of instructions videos available via the internet, e.g. under https://www.youtube.com/watch?v=sHGqwF2s-tM. For this reason, there are no further basic instruction into *Octave* in this work.