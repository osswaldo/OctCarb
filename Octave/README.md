# Octave
In this section, an instruction for instll *Octave* under Windows, macOS and Linux. In addition, a an example script for the refinement of measured WAXS/WANS data using the model of Ruland and Smarsly (2002).

All content in this repository is manly based on the following publications:
* Ruland, W. and Smarsly, B. M. (2002), X-ray scattering of non-graphitic carbon: an improved method of evaluation. J. Appl. Cryst., 35, 624-633, [doi:10.1107/S0021889802011007](https://doi.org/10.1107/S0021889802011007)
* Faber, K., Badaczewski, F., Oschatz, M., Mondin, G., Nickel, W., Kaskel, S., Smarsly, B. M. (2014), In-Depth Investigation of the Carbon Microstructure of Silicon Carbide-Derived Carbons by Wide-Angle X-ray Scattering, J. Phys. Chem. C., 118, 29, 15705-15715, [doi:10.1021/jp502832x](https://doi.org/10.1021/jp502832x)
Pfaff, T., Simmermacher, M. & Smarsly, B. M. (2018), *CarbX*: a program for the evaluation of wide-angle X-ray scattering data of non-graphitic carbons, J. Appl. Cryst., 51, 219-229, [doi:10.1107/S1600576718000195](https://doi.org/10.1107/S1600576718000195)
* Pfaff, T., Simmermacher, M. & Smarsly, B. M. (2018), *CarbX*: a program for the evaluation of wide-angle X-ray scattering data of non-graphitic carbons, J. Appl. Cryst., 51, 219-229, [doi:10.1107/S1600576718000195](https://doi.org/10.1107/S1600576718000195)

## Fit-routine.m
The file *Fit-routine.m* contains an *Octave*-script for the refinement of measured WAXS/WANS data using the model of Ruland and Smarsly (2002). The input parameters are descriped in the file itself and in https://github.com/osswaldo/ngcs/README.md.

### Output of *Fit-routine.m*
During the whole refinement, *Octave* shows a plot with the actual state of refinement (if acivated).
After every refinement step, *Octave* will return the refined values, the reason of the refinement stop, a matrix containing the refined values and the standard deviations. After the whole refinmenent, *Octave* will also return all refined parameters and the resulting microstructure parameters including the mathematically calculated erros. In addition, the used reinfment script, the plots (if activated), the refined parameters and an table containing all parameters and the refined iObs-values and a file containing the output-matrices and other useful values will be exported.

#### Console output
In addition to the informative outputs such as the start and end time, the current refinement step, etc., any errors and, after each refinement, the reason for the end of the refinement, the refined values and an associated matrix are output.

##### Reasons for the end of the refinement
* *Canceled*: Canceled by user
* *Maximum number of iterations.*: The maximum number of iterations has been arrived. You have to increase the parameter *maxIter*.
* *Change in refinement parameters too small.*: Obvious
* *Change in calculated function too small.*: Obvious, see *tolFun*

##### Refined values and matrices
The refined parameters and their values will be output after every refinement step by their name. The corresponding matrix also contains the values (1st column) and the error of it (2nd column, 1-sigma neighborhood of normal distribution). The 3rd and 4th columns contains the 2-sigma and 3-sigma neighborhoods of the normal distributions, respectiveley). An 1-sigma neighborhood of *10* means, that the error could not be calculated.

##### Common errors and how to resolve
TODO

#### File output
All files will be saved in the directory <*fitPath*>/<*filename*>

##### Plots
The plots of every refinement step and the error. For last step (refinement of all parameters together) a plot of the error (relative) and a logarthmic plot will be saved.

##### output_*.txt
This file contains all refined and constant microstructure parameters.

##### output_*.csv
In this file, all refined and constant microparameters are saved. An 1-sigma error of *10* means, that the error could not be calculated. In addition, the measured intensity (*I*), the refined intensity (*Ifit*) and the absolute (*errorAbs*) and relative )*errorRel*) error will be saved. This file can be opened directly with a table calculation program.

##### data_*.mat
Contains the following values/vactors/matrices:
* tolFun
* maxIter
* paramn
* paramFinished
* param1
* convergence1
* outp1
* result1
* stdabw1
* mat1
* param2
* convergence2
* outp2
* result2
* stdabw2
* mat2
* param3
* convergence3
* outp3
* result3
* stdabw3
* mat3
* param4
* convergence4
* outp4
* result4
* stdabw4
* mat4
* param5
* convergence5
* outp5
* result5
* stdabw5
* mat5
* rQuadratFit
* chiQuadratFit

## Installation
In principle,* Octave* is available for Microsoft Windows, Aplle macOS, GNU/Linux ad BSD systems. While a precompiled and executable* Octave* installer is available for Windows,* Octave* must be compiled on the system itself for other systems.
Independent of the system, some additional plugins for* Octave* have to be installed manually.

### Install* Octave* under Windows (tested on Windows 10 64 bit)
For Windows, an executable * .exe file can be downloaded directly: https://www.gnu.org/software/octave/download

#### Update* Octave* under Windows
Run 'Install* Octave* under Windows' again.

### Install* Octave* under macOS (tested on macOS Catalina 10.15.7)
Octave for macOS can be installed using a package manager like (Homebrew, MacPorts and Spack). In addition, a launcher app using AppleScript can also be created.
The most easiest way is to use Homebrew (https://brew.sh/), so only this way will be described.

#### Install* Octave* using Homebrew
The actual instruction is available under https://wiki.octave.org/Octave_for_macOS
0. Information: The entire process can take up to a few hours, but the actual working time in which YOU have to do something is only a few minutes.
1. Install Xcode via the App Store
2. Open a terminal and execute the following commands:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
brew update
brew install octave
```

##### Further tips
In addition,* Octave* can also be compiled against gnuplot instead of qt. This result might perform better, so it is recommended to run also the following commands:
0. Information: The entire process can take up to a few hours, but the actual working time in which YOU have to do something is only a few minutes.
1. Open the file ~/.octaverc and paste the in:
```text
setenv('GNUTERM','qt')
graphics_toolkit("gnuplot")
```
2. Open a terminal and execute the following commands:
```bash
sudo chown -R `whoami` /usr/local/share/ghostscript
brew link --overwrite ghostscript
brew install octave
# If install returns an error about not having a formular for octave, use the following command:
brew tap --repair
```

##### Update* Octave* under macOS using Homebrew
In addition,* Octave* can also be compiled against gnuplot instead of qt. This result might perform better, so it is recommended to run also the following commands:
0. Information: The entire process can take up to a few hours, but the actual working time in which YOU have to do something is only a few minutes.
1. Open a terminal and execute the following commands:
```bash
brew update && brew upgrade octave
```

### Install* Octave* under Linux
There are many different Linux distributions, so it is impossible to write one guide for all of them. Therefore, only the most relevant and widely used Linux distributions are considered. In principle,* Octave* can be installed with a precompiled package or by compiling on the operating system itself. Since the compilation is sometimes very complicated and requires a high level of prior knowledge, this procedure is NOT described here. The corresponding instructions can be found here: https://wiki.octave.org/Category:Installation

#### Debian and Debian-based (such as Ubuntu) (tested on raspbian, which is based on Debian buster)
Link: https://wiki.octave.org/Octave_for_Debian_systems
1. Open a terminal and execude the following commands:
```bash
sudo apt-get install octave
sudo apt-get install octave-io octave-statistics octave-structs octave-optim
sudo apt-get install liboctave-dev
# Optional
sudo apt-get install octave-doc octave-info octave-htmldoc
```

#### Red Hat Enterprise/CentOS
Link: https://wiki.octave.org/Octave_for_Red_Hat_Linux_systems
1. Open a terminal and execude the following commands:
```bash
sudo yum install epel-release
sudo yum install octave
yum install octave-devel
```

#### Arch Linux (not tested)
Link: https://wiki.octave.org/Octave_for_Arch_Linux
1. Open a terminal and execude the following command:
```bash
sudo pacman -S octave
```

#### Gentoo (not tested)
Link: https://wiki.octave.org/Octave_for_GNU/Linux#Gentoo
1. Open a terminal and execude the following command:
```bash
sudo emerge --ask sci-mathematics/octave
```

#### Fedora (not tested)
Link: https://wiki.octave.org/Octave_for_Red_Hat_Linux_systems
1. Open a terminal and execude the following command:
```bash
sudo dnf install octave
sudo dnf install octave-devel
```

#### openSUSE (not tested)
Link: https://wiki.octave.org/Octave_for_openSUSE
1. Open a terminal and execude the following command:
```bash
sudo zypper install octave
sudo zypper install octave-devel
# The following steps are optional:
# Install OpenBLAS https://www.openblas.net/
# Afterwards, open a terminal and execude the following commands and choose OpenBLAS:
/usr/sbin/update-alternatives --config libblas.so.3
/usr/sbin/update-alternatives --config liblapack.so.3
```

### Install* Octave* under BSD (not tested)
For FreeBSD and OpenBSD, precompiled packages* Octave* packages are available (https://wiki.octave.org/Octave_for_other_Unix_systems):

#### FreeBSD (not tested)
Open a terminal and execude the following command:
```bash
pkg_add -r octave
```

#### OpenBSD (not tested)
Open a terminal and execude the following command:
```bash
pkg_add octave
```

### *Important*: Install additional packages/plugins (optim)
In order to use the script to refine wide-angle X-ray/neutron scattering (WAXS/WANS) data of NGCs using the Octave-script in this repository, some additional plugins must be installed.
*optim* is used to perform the mathimatical refinement/minimalization operations. *optim* itself need *struct* and *statistics* as dependencies and *statistis* in turn uses *io*. If an error occurs during the installation of *optim*, it is usually because the dependencies mentioned are missing. These must then also be installed. 
Basically, either the automatic package installation program from Octave, which downloads and installs the package, or the manual package installation program, in which you have to download the packages manually, can be used for the installation * optim *.
In principle, manual installation should only be used if the automatic installation program has failed.

#### Automatic installation of *optim*
##### Installation of *optim*
The *optim* package is available under https://octave.sourceforge.io/optim/
To install *optim* either open the Octave-GUI or the Octave-CLI and execute the following command:
```octave
pkg install -forge optim
```

If the installation failed due to missing dependencies, you have to install *struct*, *statistics* and *io* **before** the installation of *optim*.
Aftwerwards, either open the Octave-GUI or the Octave-CLI and execute the following command again:
```octave
pkg install -forge optim
```

##### Installation of *struct*
To install *struct* either open the Octave-GUI or the Octave-CLI and execute the following command:
```octave
pkg install -struct optim
```

##### Installation of *statistics*
To install *statistics*, you first have to install *io*. To do this, either open the Octave-GUI or the Octave-CLI and execute the following commands:
```octave
pkg install -forge io
pkg install -forge statistics
```

### Manual installation of *optim*
If the automatic installation failed, you can also install *optim* and its dependenceis maunal.
1. Download *optim* and the needed dependecies (*structs*, *statistics* and *io*):
```links
https://octave.sourceforge.io/optim/
https://octave.sourceforge.io/struct/index.html
https://octave.sourceforge.io/statistics/index.html
https://octave.sourceforge.io/io/index.html
```

2. Install first *io*, then *structs* and *statistics* and last *optim*. To do this, either open the Octave-GUI or the Octave-CLI and execute the following commands.
 <*path_to_file*> is the path, there the previous downloaded files are saved. The quotes (') are only necessarry, if the path or file name contains spaces. <*version*> is the version of the downloaded package, i.e. the filename.
```octave
pkg install '<*path_to_file*>/io-<*version*>.tar.gz'
pkg install '<*path_to_file*>/struct-<*version*>.tar.gz'
pkg install '<*path_to_file*>/statistics-<*version*>.tar.gz'
pkg install '<*path_to_file*>/optim-<*version*>.tar.gz'
```

Copyright (C) 2021 Oliver Osswald