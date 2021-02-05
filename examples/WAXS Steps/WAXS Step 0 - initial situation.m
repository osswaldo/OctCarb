clear;

timeStart = time();

#Initial values for fitting. k & const1 should be roughly adjusted manually. nu depends on the experiment, for Cu-radiation, 4 is sufficient. Mostly, g can be 0 (switch with useGradient = true/false).
u3 = 0; #Generally not refined and not outputted
mu  = 4;
beta  = 0.5;
a3  = 3.5;
da3  = 0.4;
a3min  = a3-da3;
sig3  = 0.25;
eta  = 1;
nu  = 4;
alpha  = 0.2;
lcc  = 1.412;
sig1  = 0.1;
q  = 0;
dan  = 0;
k  = 500;
const1 = 0;
const2 = 0;

#Switch for usage of gradient g and concentrations of impurities
useGradient = false;
g  = 0;

cno = 0.0; #Concentration of disordered sp3 carbon
cH  = 0.0; #Concentration of hydrogen
cN  = 0.0; #Concentration of nitrogen
cO  = 0.0; #Concentration of oxygen
cS  = 0.0; #Concentration of sulfur

#Swtich for show a plot using the values above. Works only, if "shouldPlot = true".
plotOnly = true;

#Graphical output (has to be "false" if using octave-cli).  ('global' can be ignored, but must be present, it is necessary)
#Global variables can only be resetted restarting Octave.
global shouldPlot = true;

#Name of the series and id of the sample
name = "name";
global id = "id";
#Filename and path the currently used file, must also contain iObs.oct. The path must be changed twice.
filename = "filename.m";
#The '/' symbol must be used in the paths						 
path = '<path_to_filename>';
cd '<path_to_filename>';

#Measurement data file
measFile = '<path_to_measurement_file>';

#Corrections for Wide-Angle Neutron Scattering (WANS) experiments, only meaningful, if radiation = 1 (means neutrons scattering)
neutronCorrection = false;

#Useful for samples containing hydrogen. Using this method, a Pseudo-Voigt function will be used to fit the background instead of a polynomial. If false, the Placzek correction (a*x^2 + b) will be used for the background determination.
neutronCorrectionVoigt = false;

#Wavelength and type of radiation (0 = X-ray, 1 = neutrons)
wavelength = 1.5418;
radiation = 0;

#Calculate coherent and/or incoherent radiation. Incoherent radiation is only considered for X-ray radiation
coh = true;
inc = true;

#Type of x-values:
#theta = Scattering angle theta in °
#thetaRad = Scattering angle theta in rad.
#twoTheta = Scattering angle 2 theta in °
#twoThetaRad = Scattering angle 2 theta in rad.
#scatS = Scattering vector s = 2*sin(theta)/wavelength
#scatQ = Scattering vector q = 2*Pi*s
type = "twoTheta";

#Octave must be restarted after changing one of the below lines
#Skip n points at start
nStart = 0;

#Skip n points at the end
nEnd = 0;

#Calculate only every n point
nSkip = 1;

#Shifting the raw data up/down
nUp = 0;