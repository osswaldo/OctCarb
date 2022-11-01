# Purpose of the software
*OctCarb* is a script for ***Oct**ave* that can be used to analyze and evaluate scattering data of non-graphitic **carb**ons. Non-graphitic carbons (NGCs) are an important class of material used in acadameical research, industrial development and commercial applications. In the most of the cases, the microstructure is directly related to the macroscopical physical properties such as the high electric conductivity or chemical resitance. Using wide-angle X-ray/neutron scattering (WAXS/WANS) data or the resulting pair-distribution-function (PDF), the microstructure can be determined by using the model of Ruland & Smarsly (2002). The PDF analysis is still in progress and not maintained.

All content in this repository is manly based on the following publications:
* Ruland, W. and Smarsly, B. M. (2002), X-ray scattering of non-graphitic carbon: an improved method of evaluation. J. Appl. Cryst., 35, 624-633, [doi:10.1107/S0021889802011007      ](https://doi.org/10.1107/S0021889802011007)      
* Faber, K., Badaczewski, F., Oschatz, M., Mondin, G., Nickel, W., Kaskel, S., Smarsly, B. M. (2014), In-Depth Investigation of the Carbon Microstructure of Silicon Carbide-Derived Carbons by Wide-Angle X-ray Scattering, J. Phys. Chem. C., 118, 29, 15705-15715, [doi:10.1021/jp502832x      ](https://doi.org/10.1021/jp502832x)      
* Pfaff, T., Simmermacher, M. & Smarsly, B. M. (2018), *CarbX*: a program for the evaluation of wide-angle X-ray scattering data of non-graphitic carbons, J. Appl. Cryst., 51, 219-229, [doi:10.1107/S1600576718000195      ](https://doi.org/10.1107/S1600576718000195)      
* Osswald, O., Smarsly, B. M. (2022), J. Appl. Cryst., in preparation

# Quick start Guide:
See https://github.com/osswaldo/NGCs/blob/master/Qick-Start-Guide.md

# Where is what to find
This repository contains instruction videos, installation tutorials als examples refinement for the analysis of wide-angle X-ray/neutron (WAXS/WANS) scattering data of non-graphitic carbons (NGCs) using the Ruland & Smarsly (2002) algorithm.

## Example refinement
Hence, the refinement scripts for *Octave* as well as the corresponding "iObs.oct" file can be found. Normally, "WAXS Fit-Routine-IUCr.m" should be used and only if the error bars from the refinement parameters are from interest, but cannot be calculated, try the script "WAXS Fit-Routine-IUCr with additional step all-without-normalization.m".

### Good fit
An example for a "good" *Octave* fit with the important fit regions marked in green (important), yellow (nice to have) and red (unimportant).

### WAXS Steps
An example refinement with given scattering data and examples/scripts for each refinement step.

## Instruction Videos
As the name says, some instruction videos for the installation of *Octave*, fitting of WAXS/WANS data of NGCs and the compilation of an iObs.oct-file.

## *Octave*
The directory *Octave* contains instructions for installing *GNU Octave* on Windows, macOS and Linux. In addition, a refinment script for get the microstructure date of non-graphitic carbons (NGCs) from wide-angle X-ray/neutron (WAXS/WANS) scattering data using the Ruland & Smarsly (2002) algorithm. Also, an instruction for compiling the needed *iObs.oct* file as well as some already compiled *iObs.oct* files are also provided, which can be used to call *iObs* from Octave. For optimal functionality, the file should be renamed to *iObs.oct* after downloading/compiling.
If you never worked with* Octave* so far, it is highly recommended to take a look on an instruction video. Due to the high distribution and availability of *Octave*, there are a lot of instructions videos available via the internet, e.g. under https://www.youtube.com/watch?v=sHGqwF2s-tM. For this reason, there are no further basic instruction into *Octave* in this work.

##Useful-links
See https://github.com/osswaldo/NGCs/tree/master/Useful-links.md

Copyright (C) 2022 Oliver Osswald