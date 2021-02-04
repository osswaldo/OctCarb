# NGCs
Non-graphitic carbons (NGCs) are an important class of material used in acadameical research, industrial development and commercial applications. In the most of the cases, the microstructure is directly related to the macroscopical physical properties such as the high electric conductivity or chemical resitance. Using wide-angle X-ray/neutron scattering (WAXS/WANS) or the resulting pair-distribution-function (PDF), the microstructure can be determined by using the model of Ruland & Smarsly (2002). The PDF analysis is still in progress and can be used as soon as possible.

All content in this repository is manly based on the following publications:
* Ruland, W. and Smarsly, B. M. (2002), X-ray scattering of non-graphitic carbon: an improved method of evaluation. J. Appl. Cryst., 35, 624-633, [doi:10.1107/S0021889802011007](https://doi.org/10.1107/S0021889802011007)
* Faber, K., Badaczewski, F., Oschatz, M., Mondin, G., Nickel, W., Kaskel, S., Smarsly, B. M. (2014), In-Depth Investigation of the Carbon Microstructure of Silicon Carbide-Derived Carbons by Wide-Angle X-ray Scattering, J. Phys. Chem. C., 118, 29, 15705-15715, [doi:10.1021/jp502832x](https://doi.org/10.1021/jp502832x)
Pfaff, T., Simmermacher, M. & Smarsly, B. M. (2018), *CarbX*: a program for the evaluation of wide-angle X-ray scattering data of non-graphitic carbons, J. Appl. Cryst., 51, 219-229, [doi:10.1107/S1600576718000195](https://doi.org/10.1107/S1600576718000195)
* Pfaff, T., Simmermacher, M. & Smarsly, B. M. (2018), *CarbX*: a program for the evaluation of wide-angle X-ray scattering data of non-graphitic carbons, J. Appl. Cryst., 51, 219-229, [doi:10.1107/S1600576718000195](https://doi.org/10.1107/S1600576718000195)

## iObs
The directory *iObs* contains a C++ code for calculating theroetical intensity of wide-angle X-ray/neutron (WAXS/WANS) scattering of non-graphitic carbons (NGCs) using the Ruland & Smarsly (2002) algorithm. In addition, the pair-distribution-function (PDF) can also be calculated either using WAXS/WANS data and a fourier transformation ("pdf") or calculating iObs and PDF and one step ("iObsPDF").

## Octave
The directory *Octave* contains instructions for installing *GNU Octave* on Windows, macOS and Linux. In addition, a refinment script for get the microstructure date of non-graphitic carbons (NGCs) from wide-angle X-ray/neutron (WAXS/WANS) scattering data using the Ruland & Smarsly (2002) algorithm. Also, an instruction for compiling the needed *iObs.oct* file as well as some already compiled *iObs.oct* files are also provided, which can be used to call *iObs* from Octave. For optimal functionality, the file should be renamed to *iObs.oct* after downloading/compiling.
If you never worked with* Octave* so far, it is highly recommended to take a look on an instruction video. Due to the high distribution and availability of *Octave*, there are a lot of instructions videos available via the internet, e.g. under https://www.youtube.com/watch?v=sHGqwF2s-tM. For this reason, there are no further basic instruction into *Octave* in this work.

## In- and Output parameters
In general, the input parameters for *iObs*, *iObsPDF* and the associated *.oct* files are the same. The only difference between *iObs* and *iObsPDF* is that *iObsPDF* requires some additional parameters.
Further information about all parameters can be found in the pulblications on top.

### Input parameters
#### iObs
| Number 	| Parameter            	| Description                                                                                 	| Type   	| Default value/example for *Octave*                            	|
|--------	|----------------------	|---------------------------------------------------------------------------------------------	|--------	|---------------------------------------------------------------	|
| 1      	| cno                  	| Concentration   of unorganized carbon                                                       	| double 	| cno = 0;                                                      	|
| 2      	| mu                   	| Shape factor for gamma function to calculate N, Lc and kapc                                 	| double 	| mu = 4;                                                       	|
| 3      	| beta                 	| Shape factor for gamma function to calculate N, Lc and kapc                                 	| double 	| beta = 0.5;                                                   	|
| 4      	| a3                   	| Average layer distance                                                                      	| double 	| a3 = 0.5;                                                     	|
| 5      	| da3                  	| Minimal layer distance                                                                      	| double 	| da3 = 0.4; a3min = a3-da3;                                    	|
| 6      	| sig3                 	| Standard deviation of a3                                                                    	| double 	| sig3 = 0.25;                                                  	|
| 7      	| u3                   	| Thermal motion                                                                              	| double 	| u3 = 0;                                                       	|
| 8      	| eta                  	| Homogeneity of the stacks                                                                   	| double 	| eta = 1;                                                      	|
| 9      	| nu                   	| Shape factor for gamma function to calculate La, Lm and kapa                                	| int    	| nu = 4; (for Cu-radiation) nu = 7; (for neutron-radiaiton)    	|
| 10     	| alpha                	| Shape factor for gamma function to calculate La, Lm and kapa                                	| double 	| alpha = 0.2;                                                  	|
| 11     	| lcc                  	| Average C-C bond length                                                                     	| double 	| lcc = 1.412;                                                  	|
| 12     	| sig1                 	| Standard deviation of lattice constant a = √3 \* lcc                                        	| double 	| sig1 = 0.1;                                                   	|
| 13     	| q                    	| Preferred orientation                                                                       	| double 	| q = 0;                                                        	|
| 14     	| cH                   	| Concentration of unorganized hydrogen                                                       	| double 	| cH = 0;                                                       	|
| 15     	| cN                   	| Concentration of unorganized nitrogen                                                       	| double 	| cN = 0;                                                       	|
| 16     	| cO                   	| Concentration of unorganized oxygen                                                         	| double 	| cO = 0;                                                       	|
| 17     	| cS                   	| Concentration of unorganized sulfur                                                         	| double 	| cS = 0;                                                       	|
| 18     	| dan                  	| Anisotropy of atomic form factor of carbon                                                  	| double 	| dan = 0;                                                      	|
| 19     	| k                    	| Normalization constant for log10(k \* Ie.u. + const1) + const2                              	| double 	| k = 500; (depending on radiation intensity and sample amount) 	|
| 20     	| const1               	| Constant shift for log10(k \* Ie.u. + const1) + const2                                      	| double 	| const1 = 0;                                                   	|
| 21     	| const2               	| Non-constant (linear) shift for log10(k \* Ie.u. + const1) + const2                         	| double 	| const2 = 0;                                                   	|
| 22     	| useQ                 	| Switch for additional parameter for incoherent background                                   	| bool   	| useQ = false;                                                 	|
| 23     	| b                    	| Additional parameter for incoherent background                                              	| double 	| b = 0;                                                        	|
| 24     	| useA                 	| Switch for additional adsorption parameters                                                 	| bool   	| useQ = true;                                                  	|
| 25     	| density              	| Density of the sample in g/cm^3                                                             	| double 	| density = 2.2;                                                	|
| 26     	| sampleThickness      	| Thickness of the sample in cm                                                               	| double 	| sampleThickness = 0.3;                                        	|
| 27     	| transmission         	| Transmisison geometry (if false, reflection geometry is assumed)                            	| bool   	| transmission = false;                                            	|
| 28     	| absorptionCorrection 	| Additional correction factor for the calculation absorption coefficient                     	| double 	| absorptionCorrection = 1;                                     	|
| 29     	| useP                 	| Switch for addiational polarization correction                                              	| bool   	| useP = true;                                                  	|
| 30     	| polarizedBeam        	| Is the beam polarized?                                                                      	| bool   	| polarizedBeam = false;                                        	|
| 31     	| polarizationDegree   	| Polarization direction of beam in °			                                              	| double 	| polarizationDegree = 0;                                       	|
| 32     	| useGradient          	| Switch for additional exponential damping                                                   	| bool   	| useQ = false;                                                 	|
| 33     	| g                    	| Factor for exponential damping of the scattering intensity with Ie.u.’ = exp(g\*s) \* Ie.u. 	| double 	| g = 0;                                                        	|
| 34     	| useCorrAutoColl      	| Switch for additional conversion from fixed irradiated length to fixed slit                 	| bool   	| useQ = false;                                                 	|
| 35     	| par_r                	| Radius of the goniometer (in cm; fixed due to experiment)                                   	| double 	| par_r = 14; (depends on the experiment)                       	|
| 36     	| par_delta            	| Divergance angle (in °; to choose by user)                                                  	| double 	| par_delta = 4;                                                	|
| 37     	| par_l                	| Irradiated length (in cm; fixed during measurement)                                         	| double 	| par_l = 7; (depends on the experiment)                        	|
| 38     	| radiationType        	| Type of radiation (0 = X-ray, 1 = neutron)                                                  	| int    	| radiation = 0;                                                	|
| 39     	| wavelength           	| Wavelength                                                                                  	| double 	| wavelength = 1.5418; (Cu-radiation)                           	|
| 40     	| S                    	| Vector, which contains the scattering vector values (s =  sin(θ) / λ)                       	| vector 	| s = [0; 0.01; 0.02; ...];                                     	|
| 41     	| coh                  	| Switch for calculating coherent scattering                                                  	| bool   	| coh = true;                                                   	|
| 42     	| inc                  	| Switch for calculating incoherent scattering                                                	| bool   	| inc = true; (only meaningful for X-ray radiation)             	|


#### iObsPDF
| Number 	| Parameter            	| Description                                                                                               	| Type   	| Default value/example for *Octave*                            	|
|--------	|----------------------	|-----------------------------------------------------------------------------------------------------------	|--------	|---------------------------------------------------------------	|
| 1      	| cno                  	| Concentration   of unorganized carbon                                                                     	| double 	| cno = 0;                                                      	|
| 2      	| mu                   	| Shape factor for gamma function to calculate N, Lc and kapc                                               	| double 	| mu = 4;                                                       	|
| 3      	| beta                 	| Shape factor for gamma function to calculate N, Lc and kapc                                               	| double 	| beta = 0.5;                                                   	|
| 4      	| a3                   	| Average layer distance                                                                                    	| double 	| a3 = 0.5;                                                     	|
| 5      	| da3                  	| Minimal layer distance                                                                                    	| double 	| da3 = 0.4; a3min = a3-da3;                                    	|
| 6      	| sig3                 	| Standard deviation of a3                                                                                  	| double 	| sig3 = 0.25;                                                  	|
| 7      	| u3                   	| Thermal motion                                                                                            	| double 	| u3 = 0;                                                       	|
| 8      	| eta                  	| Homogeneity of the stacks                                                                                 	| double 	| eta = 1;                                                      	|
| 9      	| nu                   	| Shape factor for gamma function to calculate La, Lm and kapa                                              	| int    	| nu = 4; (for Cu-radiation) nu = 7; (for neutron-radiaiton)    	|
| 10     	| alpha                	| Shape factor for gamma function to calculate La, Lm and kapa                                              	| double 	| alpha = 0.2;                                                  	|
| 11     	| lcc                  	| Average C-C bond length                                                                                   	| double 	| lcc = 1.412;                                                  	|
| 12     	| sig1                 	| Standard deviation of lattice constant a = √3 \* lcc                                                      	| double 	| sig1 = 0.1;                                                   	|
| 13     	| q                    	| Preferred orientation                                                                                     	| double 	| q = 0;                                                        	|
| 14     	| cH                   	| Concentration of unorganized hydrogen                                                                     	| double 	| cH = 0;                                                       	|
| 15     	| cN                   	| Concentration of unorganized nitrogen                                                                     	| double 	| cN = 0;                                                       	|
| 16     	| cO                   	| Concentration of unorganized oxygen                                                                       	| double 	| cO = 0;                                                       	|
| 17     	| cS                   	| Concentration of unorganized sulfur                                                                       	| double 	| cS = 0;                                                       	|
| 18     	| dan                  	| Anisotropy of atomic form factor of carbon                                                                	| double 	| dan = 0;                                                      	|
| 19     	| k                    	| Normalization constant for log10(k \* Ie.u. + const1) + const2                                            	| double 	| k = 500; (depending on radiation intensity and sample amount) 	|
| 20     	| const1               	| Constant shift for log10(k \* Ie.u. + const1) + const2                                                    	| double 	| const1 = 0;                                                   	|
| 21     	| const2               	| Non-constant (linear) shift for log10(k \* Ie.u. + const1) + const2                                       	| double 	| const2 = 0;                                                   	|
| 22     	| useQ                 	| Switch for additional parameter for incoherent background                                                 	| bool   	| useQ = false;                                                 	|
| 23     	| b                    	| Additional parameter for incoherent background                                                            	| double 	| b = 0;                                                        	|
| 24     	| useA                 	| Switch for additional adsorption parameters                                                               	| bool   	| useQ = true;                                                  	|
| 25     	| density              	| Density of the sample in g/cm^3                                                                           	| double 	| density = 2.2;                                                	|
| 26     	| sampleThickness      	| Thickness of the sample in cm                                                                             	| double 	| sampleThickness = 0.3;                                        	|
| 27     	| transmission         	| Transmisison geometry (if false, reflection geometry is assumed)                                          	| bool   	| transmission = false;                                                 	|
| 28     	| absorptionCorrection 	| Additional correction factor for the calculation absorption coefficient                                   	| double 	| absorptionCorrection = 1;                                     	|
| 29     	| useP                 	| Switch for addiational polarization correction                                                            	| bool   	| useP = true;                                                  	|
| 30     	| polarizedBeam        	| Is the beam polarized?                                                                                    	| bool   	| polarizedBeam = false;                                        	|
| 31     	| polarizationDegree   	| Polarization direction of beam in °				                                                           	| double 	| polarizationDegree = 0;                                       	|
| 32     	| useGradient          	| Switch for additional exponential damping                                                                 	| bool   	| useQ = false;                                                 	|
| 33     	| g                    	| Factor for exponential damping of the scattering intensity with Ie.u.’ = exp(g\*s) \* Ie.u.               	| double 	| g = 0;                                                        	|
| 34     	| useCorrAutoColl      	| Switch for additional conversion from fixed irradiated length to fixed slit                               	| bool   	| useQ = false;                                                 	|
| 35     	| par_r                	| Radius of the goniometer (in cm; fixed due to experiment)                                                 	| double 	| par_r = 14; (depends on the experiment)                       	|
| 36     	| par_delta            	| Divergance angle (in °; to choose by user)                                                                	| double 	| par_delta = 4;                                                	|
| 37     	| par_l                	| Irradiated length (in cm; fixed during measurement)                                                       	| double 	| par_l = 7; (depends on the experiment)                        	|
| 38     	| radiationType        	| Type of radiation (0 = X-ray, 1 = neutron)                                                                	| int    	| radiation = 0;                                                	|
| 39     	| wavelength           	| Wavelength                                                                                                	| double 	| wavelength = 1.5418; (Cu-radiation)                           	|
| 40     	| S                    	| Vector, which contains the scattering vector values (s =  sin(θ) / λ)                                     	| vector 	| s = [0; 0.01; 0.02; ...];                                     	|
| 41     	| coh                  	| Switch for calculating coherent scattering                                                                	| bool   	| coh = true;                                                   	|
| 42     	| inc                  	| Switch for calculating incoherent scattering                                                              	| bool   	| inc = true; (only meaningful for X-ray radiation)             	|
| 43     	| Smin                 	| Minimum scattering vector for calculating *iObs*, which is required for calculating the PDF values        	| double 	| Smin = 0.1;                                                   	|
| 44     	| Smax                 	| Maximal scattering vector for calculating *iObs*, which is required for calculating the PDF values        	| double 	| Smax = 25;                                                    	|
| 45     	| Sstep                	| Step-width for scattering vector for calculating *iObs*, which is required for calculating the PDF values 	| double 	| Sstep = 0.01;                                                 	|
| 46     	| Rmin                 	| Minimum atomic distance for the PDF-values                                                                	| double 	| Rmin = 1;                                                     	|
| 47     	| Rmax                 	| Maximum atomic distance for the PDF-values                                                                	| double 	| Rmax = 40;                                                    	|
| 48     	| Rstep                	| Step-width for atomic distance for the PDF-values                                                         	| double 	| Rstep = 0.1;                                                  	|
| 49     	| dfactor              	| D-factor for calculating the PDF-values (identical influence like *scale*)                                	| double 	| dfactor = 1;                                                  	|
| 50     	| scale                	| Scale-factor for calculating the PDF-values (identical influence like *dfactor*)                          	| double 	| scale = 1;                                                    	|


#### Relationship to microstructure parameters


### Output parameters
Both, *iObs* and *iObsPDF* return a matrix that contains the x/y values of the calculations.

#### iObs
* First column: Scattering vector values (*s* in Å^-1)
* Second column: Oberseved Intensity (*iObs*)

#### iObsPDF
* First column: Atomic distance (*r* in Å)
* Second column: Local electron density (*G(r)*)

## Common warnings/errors
There are several errors, which can occur during the refinement process. The most common and their solutions are written down here:

### warning: lower and upper bounds identical for some parameters, fixing the respective parameters
As the warning says, for one or more parameters, the upper and lower bound are identical. You can check the parameters and ignore the warning.

### warning: matrix singular to machine precision, rcond = nan
Some values in the matrix are too close to 0 and your PC cannot calculate it. It is not a big problem, the PC will use 0 instead and forward the refinement. In general, this warning can be ignored.

### some lower bounds larger than upper bounds
You set some lower bound to a higher value than the upper bound. You have to check and correct the parameters.

### no free parameters
The lower and upper bounds of all parameters are fixed to a certain value, so there is nothing which can be refined. At least one parameter must be free.

### __plt2vv__: vector lengths must match
This error typical occurs, if a parameter is changed, which influences the raw data, e.g., measFile, type, nStart, nEnd, nSkip. In the most of the cases, the error can resolved by restarting Octave.

### svd: cannot take SVD of matrix containing Inf or NaN values
There are a lot of reasons, why this error is shown. It might be, that the starting parameters are too “bad”, the amount of measurement is to low, i.e. nSkip is to high or there is to less RAM available. First, you can ignore this error,* Octave* will try to calculate the Matrix again (10 times). If the script stops, the start parameters should be refined manually.

### error: Too many errors when refine
This is not a specific error, it occurs, if something goes wrong. There should be another error or warning in the output and this is the reason, why the script stops. In principle, all manually changed parameters should be checked.

### The values of the error Matrix are 10, 20 and 30
These values are not the real errors, but something like a placeholder or default value. A standard deviation of 10 means, that the error cannot be calculated correctly. This is a very common mistake when trying to refine too many parameters at the same time.

Copyright (C) 2021 Oliver Osswald