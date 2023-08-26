# Compile und usage of an *iObs.oct*-file
To refine wide-angle X-ray/neutron scattering (WAXS/WANS) of non-graphitic carbons (NGCs) using Octave, the theoretical observed intensity (iObs) has to be calculated in Octave. This is done with the C ++ *iObs* code, which has to be compiled into an *iObs.oct* file.
However, it is also possible to us an *iObsPDF.oct* file to calculate pair-distribution-function (PDF) data instead of WAXS/WANS data. A high maximum scattering vector (*s* >= 20 Å^-1), nu >= 7 and the usage of neutron radiation are recommended for this.
The files can be either downloaded from this page or compiled by yourself using the follwing instructions.
Please note that the files already compiled here may not always be up-to-date and will not work on all systems.

## Example files
The provided example oct files are compiled with Octave 5.2.0 (Windows versions also 6.3.0). They may also work with other Octave/operating system versions, but this has not been tested.

## Numerical recipes
iObs uses parts of the Numerical Recipes (written in C), available under a commercial license. The machine-readable code of CarbX must therefore not be distributed. However, some public domain Numerical Recipes code is included in source code archive.

### Included
* nrcomplex.h
* nrutil.h
* nrcomplex.c
* nrutil.c

### Missing for iObs
* bessj0.c
* bsstep.c
* factln.c
* fpoly.c
* gammln.c
* gaulag.c
* hypdrv.c
* hypgeo.c
* hypser.c
* midinf.c
* mmid.c
* nr.h
* odeint.c
* pzextr.c
* qsimp.c
* spline.c
* splint.c
* trapzd.c

### Missing for iObsPDF (additional)
* fleg.c
* pythag.c
* ran1.c
* svbksb.c
* svdcmp.c
* svdfit.c
* svdvar.c

## Compile
To compile an *.oct* file, *Octave* must be installed first and the *mkoctfile* command must be executed second. The following commands are only examples and in most cases the paths need to be adjusted.

### Microsoft Windows
1. Install* Octave* as described in https://github.com/osswaldo/ngcs/Octave/README.md
2. To compile *iObs* to *iObs.oct*, just open the command line and execute the following commands (the paths may need to be adjusted):
```cmd
cd C:\Octave\Octave-6.3.0\mingw64\bin
mkoctfile -L..\lib\octave\6.3.0 -I..\include\octave-6.3.0\octave "<*path_to_iObs*>\iObs.cpp"
mkoctfile -L..\lib\octave\6.3.0 -I..\include\octave-6.3.0\octave "<*path_to_iObsPDF*>\iObsPDF.cpp"
```
This generates an *iObs.oct* and an *iObsPDF.oct* file in the directory *C:\Octave\Octave-6.3.0\mingw64\bin\mkoctfile*.

### Apple macOS
1. Install Octave
To compile *iObs* to *iObs.oct*, just open the command line and execute the following commands (the paths may need to be adjusted, <*path_to_iObs*> must be replaced):
```bash
mkoctfile -I/usr/local/bin/octave "<*path_to_iObs*>/iObs.cpp"
mkoctfile -I/usr/local/bin/octave "<*path_to_iObsPDF*>/iObsPDF.cpp"
```
This generates an *iObs.oct* and an *iObsPDF.oct* file in the current directory.

### Linux
1. Install Octave
To compile *iObs* to *iObs.oct*, just open the command line and execute the following commands (the paths may need to be adjusted, <*path_to_iObs*> must be replaced):
```bash
mkoctfile -I/usr/local/bin/octave "<*path_to_iObs*>/iObs.cpp"
mkoctfile -I/usr/local/bin/octave "<*path_to_iObsPDF*>/iObsPDF.cpp"
```
This generates an *iObs.oct* and an *iObsPDF.oct* file in the current directory.

Copyright (C) 2023 Oliver Osswald