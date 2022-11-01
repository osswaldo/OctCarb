# Quick start guide
## Key points (ALWAYS do)
1. Download Windows-Installation.zip move/copy the "NGC Analysis" folder to the drive under "C:\". Please pay attention to 32/64 bit version. For reading, see here: https://support.microsoft.com/de-de/windows/32-bit-und-64-bit-windows-h%C3%A4ufig-posed-questions-c6ca9541-8dce-4d48-0415 -94a3faa2e13d

2. Install Octave (only versions 6.3.0 and 5.2.0 was tested)

3. Open and run the Setup.m script. Only the "optim" package is installed and tested here using the "pkg install -forge optim" command.

4. Copy the "Fit-Routine.m" file and rename it by sample name/ID (to avoid overwriting old data)

5. Enter a name (=sample series) and an ID (=exact sample) in lines 53 and 54. The results are then stored in the "name/id" subfolder.

6. Only if a path other than "C:/NGC-Analysis" is used: change this in lines 58 and 59 (pay attention to "/" as folder separator).

7. Enter the path to the measurement data in line 62 (important: do NOT use a header, use data directly). By default, 2Theta is selected as the X-axis and a linear Y-axis. Other X-axes can be selected in line 85.

8. VERY IMPORTANT (lines 142 - 145): (Slit correction for from a variable slit (=fixed irradiated length) to a fixed slit (=variable irradiated length) (default: fixed slit)))!!!

9. Click on "Save file and run" (yellow arrow in gray gear wheel) at the top

10. Wait and view results in the "Command Window" or under "C:\NGC-Analysis\name\id".



## Miscellaneous:
11. If something is known about the sample (e.g. from results at a similar temperature), the starting values ​​in lines 18-33 can (but don't have to!) be changed.

12. Additional exponential damping of the data. Should not be used("= false"), only in special cases.

13. If known, give the results of the elemental analysis or the fraction of carbon sp3 (e.g. from XPS) in lines 39-43.

14. If no fit is to be carried out, but only a sample plot is to be simulated: In line 46 "plotOnly = true;" put.

15. Enter the path to the measurement data in line 62 (important: do NOT use a header, use data directly). By default, 2Theta is selected as the X-axis and a linear Y-axis. Other X-axes can be selected in line 85.

16. If neutron scattering is used: In lines 65 and 68 indicate whether the background should be taken into account automatically and whether a Voigt function (with hydrogen content > 0.5 wt%) should be used for this.

17. If radiation other than copper-K alpha (1.5418 A) is used: State this on line 71.

18. If neutron scattering is used: In line 72 "radiation = 1;" specify.

19. In lines 75 and 76 only coherent and incoherent scattering can be selected. This can be used to view the amount of background scattering or to fit it manually (e.g. when there is a high proportion of foreign atoms).

20. In lines 89 - 98 you can select whether points at the beginning or end should be left out or the measuring points should be skipped (= skipped) in order to speed up the fit.

21. In lines 102 - 115 the lower/upper limits for the parameters to be adjusted can be selected. Normally you don't need to change anything here.

22. In lines 118 - 127 the maximum number of fit iterations and the fit accuracy as well as the weighting of the data points can be selected. You may have to change it for other devices with a larger spread or more/fewer points, only "trial and error" helps here.

23. In lines 129 - 140, a Debye-Waller factor (default: off, since the effect is too small), an absorption correction (default: 3mm sample thickness and transmission geometry) and a beam polarization (default: unpolarized) can be made.