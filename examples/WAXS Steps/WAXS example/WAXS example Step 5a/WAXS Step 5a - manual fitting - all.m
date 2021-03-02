clear;

timeStart = time();

#Initial values for fitting. k & const1 should be roughly adjusted manually. nu depends on the experiment, for Cu-radiation, 4 is sufficient. Mostly, g can be 0 (switch with useGradient = true/false).
u3 = 0; #Generally not refined and not outputted
mu     = 8;        # Step 3a
beta   = 3;        # Step 3a
a3     = 3.5;      # Step 3a
da3    = 0.6;      # Step 3a
a3min  = a3-da3;
sig3   = 0.3;      # Step 3a
eta    = 0.87;     # Step 3a
nu     = 4;        # Step 1
alpha  = 0.3;      # Step 4a
lcc    = 1.412;    # Step 4a
sig1   = 0.15;     # Step 4a
q      = 0;        # Step 3a
dan    = 0;
k      = 780;      # Step 2     # Step 3a     # Step 4a
const1 = 1500;     # Step 2     # Step 3a     # Step 4a
const2 = 0;

#Switch for usage of gradient g and concentrations of impurities
useGradient = false;
g      = 0;

cno = 0.0; #Concentration of disordered sp3 carbon
cH  = 0.0; #Concentration of hydrogen    # Step 2 (known from elemental analysis)
cN  = 0.0; #Concentration of nitrogen    # Step 2 (known from elemental analysis)
cO  = 0.05; #Concentration of oxygen     # Step 2 (known from elemental analysis)
cS  = 0.02; #Concentration of sulfur     # Step 2 (known from elemental analysis)

#Swtich for show a plot using the values above. Works only, if "shouldPlot = true".
plotOnly = true;

#Graphical output (has to be "false" if using octave-cli).  ('global' can be ignored, but must be present, it is necessary)
#Global variables can only be resetted restarting Octave.
global shouldPlot = true;

#Name of the series and id of the sample
name = "WAXS example";                             # Step 1
global id = "WAXS example Step 5a";                # All steps
#Filename and path the currently used file, must also contain iObs.oct. The path must be changed twice.
#The '/' symbol must be used in the paths
filename = "WAXS Step 5a - manual fitting - all.m";                                             # All steps
path = 'D:/OneDrive/Github/NGCs/examples/WAXS Steps';                      # Step 1
cd 'D:/OneDrive/Github/NGCs/examples/WAXS Steps';                          # Step 1

#Measurement data file
measFile = 'D:/OneDrive/Github/NGCs/examples/WAXS Steps/WAXS example data.xy'; # Step 1

#Corrections for Wide-Angle Neutron Scattering (WANS) experiments, only meaningful, if radiation = 1 (means neutrons scattering)
neutronCorrection = false;

#Useful for samples containing hydrogen. Using this method, a Pseudo-Voigt function will be used to fit the background instead of a polynomial. If false, the Placzek correction (a*x^2 + b) will be used for the background determination.
neutronCorrectionVoigt = false;

#Wavelength and type of radiation (0 = X-ray, 1 = neutrons)
wavelength = 1.5418;                               # Step 1
radiation = 0;                                     # Step 1

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
type = "twoTheta";                                 # Step 1

#Octave must be restarted after changing one of the below lines
#Skip n points at start
nStart = 0;

#Skip n points at the end
nEnd = 0;

#Calculate only every n point
nSkip = 1;     # Step 1		# Step 5a

#Shifting the raw data up/down
nUp = 0;

#Lower and upper bounds for fitting parameters (can be changed, but does not have to be)
#mu, beta, a3, da3, sig3, eta
lb1 = [0.01; 0.01; 2.50; 0; 0; 0];
ub1 = [10;   20;   5.00; 1; 1; 1];

#nu, alpha, sig1
lb2 = [nu; 0.01; 0];
ub2 = [nu; 2;    2];

#q, dan, k, const1, const2, g
lb3 = [q; dan; 0.0001; -100000; const2; -1];
ub3 = [q; dan;  10000; +100000; const2; +1];

#lcc
lb4 = [1.2];
ub4 = [1.8];

#Function tolerance
tolFun = 1e-10;

#Maximal iterations per fit step
maxIter = 50;

#Weight of measurement points
#normal = 1 (every point has the same weight)
#weight = 1/y (default; lower intensity has lower weight)

weight = "weight";

useQ    = false; #Additional Debye-Waller-factor
b       = 0.002; #Factor

useA                  = true;    #Absorption correction
density               = 2.2;     #Density of sample in g/cm^3
sampleThickness       = 0.3;     #Sample thickness in cm
transmission          = false;   #Transmission geometry (if false, reflection geometry is assumed)
absorptionCorrection  = 1;       #Correction factor for absorption coefficient (multiplicated to theoretical coefficient)

useP                  = true;    #Polarization correction
polarizedBeam         = false;   #Do you use a polarized beam?
polarizationDegree    = 0;       #Polarization direction of beam in °.

useCorrAutoColl = false;  #Slit correction
par_r           = 14;	    #Radius of the diffractometer (Debye-Scherrer) in cm
par_delta       = 4;	    #Divergence angle in ° (it is converted as if this fixed slit were inside)
par_l           = 5;	    #Irradiated length in cm

#Function to call von iObs.oct. Additional refinement parameters/correction terms can be activated.
function fun = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc)
  iObsOut = iObs(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc);
  fun = iObsOut(2,:);
endfunction

#Lines below should be left as it is
# # # # #
 # # # # 
# # # # #
 # # # # 
# # # # #
#Lines below should be left as it is

#Load optim package
pkg load optim;

#Vector containing start parameters
paramn = [cno; mu; beta; a3; da3; sig3; u3; eta; nu; alpha; sig1; lcc; q; dan; k; const1; const2; g];
data = load(measFile);
x = data(:,1)';
yn = data(:,2)';

switch (type)
  case "theta"
	  x = 2/wavelength * sin(x*pi/180);
  case "thetaRad"
    x = 2/wavelength * sin(x);
  case "twoTheta"
    x = 2/1.5418 * sin(x/2*pi/180);
  case "twoTheta"
    x = 2/1.5418 * sin(x/2*pi/180);
  case "twoThetaRad"
    x = 2/1.5418 * sin(x/2);
  case "scatS"
    x = x;
  case "scatQ"
    x = x/(2*pi);
  otherwise
    error("Unknown x-value type")
endswitch

for i=1:length(x)-nStart
	xNew1(i) = x(i+nStart);
	ynNew1(i) = yn(i+nStart);
endfor

x = xNew1;
yn = ynNew1;

for i=1:length(x)-nEnd
	xNew2(i) = x(i);
	ynNew2(i) = yn(i);
endfor

x = xNew2;
yn = ynNew2;
for i=1:length(x)/nSkip
	xNew3(i) = x(i*nSkip);
	ynNew3(i) = yn(i*nSkip);
endfor

x = xNew3;
yn = ynNew3;
for i=1:length(x)
	xNew4(i) = x(i);
	ynNew4(i) = yn(i)+nUp;
endfor

x = xNew4;
yn = ynNew4;

function gauss = gauss(x, k, w, x0, y0)
  for i=1:1:length(x)
    gauss(i)=y0 + k*exp(-log(2)*((x(i)-x0)/w)^2);
  endfor
endfunction

function lorentz = lorentz(x, k, w, x0, y0)
  for i=1:1:length(x)
    lorentz(i)=y0 + k/(1+((x(i)-x0)/w)^2);
  endfor
endfunction

function pseudoVoigt = pseudoVoigt(x, k, w, x0, y0, eta)
  for i=1:1:length(x)
    pseudoVoigt(i) = (1-eta)*gauss(x(i), k, w, x0, y0) + eta*lorentz(x(i), k, w, x0, y0);
  endfor
endfunction

function SofQ = placzekCorrection(x, y, minX = 0, maxX = 100)
  for i=1:length(x)
    if (minX <= x <= maxX)
      xnew(i) = x(i);
      ynew(i) = y(i);
    endif
  endfor

  M = 2;
  N = length(xnew);
  K = 2;
  EPSH = 1e-100000;
  MAXIT = 50;
  
  xnew = xnew;  
  ynew = ynew;
  
  A = polyfit(xnew, ynew, M);
  
  fy = polyval(A, xnew);
  
  fit = y - fy;
  
  minimalValue = min(fit);
  
  #Calculation SofQ and adding the minimal value of "SofQ" to prevent negative values
  SofQ = y - fy + 1 - minimalValue + 10;
endfunction

function SofQ = lorgaunCorrection(x, y, minX = 0, maxX = 100)
  for i=1:length(x)
    if (minX <= x <= maxX)
      xnew(i) = x(i);
      ynew(i) = y(i);
    endif
  endfor
  
  #options.AutoScaling = autoscaling;
  #options.FunValCheck = funValCheck;
  options.MaxIter = 50;
  options.TolFun = 1e-100000;
  
  #Normalization k, width w, ration eta
  options.lbound = [0.01; 0.1; 0];
  options.ubound = [1e10; 100; 1];
  
  pin = [50000; 15; 0.5];
  
  #Normalization k, width w, middle, x0, shift y0, ration eta
  f = @ (p, x) (pseudoVoigt(x, p(1), p(2), 0, 0, p(3)));
  
  [p, fy, cvg, outp] = nonlin_curvefit(f, pin, xnew, ynew, options);
  
  
  fy = fy;
  p = p
  cvg = cvg
  
  fit = y - fy;
  
  minimalValue = min(fit);
  
  #Calculation SofQ and adding the minimal value of "SofQ" to prevent negative values
  SofQ = y - fy + 1 - minimalValue + 10;
endfunction

if neutronCorrection == true
	if neutronCorrectionVoigt == true
    yn = lorgaunCorrection(x, yn);
  else
    yn = placzekCorrection(x, yn);
  endif
endif

s = x;
global ynglobal = yn;

wtNormal = ones(size(yn));
for i=1:length(yn)
	if weight == "normal"
		wtWeight(i) = wtNormal;
	else
		wtWeight(i) = 1/(yn(i)+10);
	endif
endfor

#OutputPath
fitPath = strcat(path, "/", name, '/', id);

#Make dir's
mkdir(path, name);
mkdir(strcat(path, "/", name), id);

#Vector containing start-parameters
paramStart = [mu; beta; a3; da3; sig3; eta; nu; alpha; sig1; lcc; q; dan; k; const1; const2; g];

#Parameters for fit steps
paramn1 = [mu; beta; a3; da3; sig3; eta];
paramn2 = [nu; alpha; sig1];
paramn3 = [q; dan; k; const1; const2; g];
paramn4 = [lcc];

#Bounds for last step (order!)
lb5 = [lb1; lb2(1); lb2(2); lb4; lb2(3); lb3];
ub5 = [ub1; ub2(1); ub2(2); ub4; ub2(3); ub3];

#Typical x-value (scattering vector s)
typicalX = 1;

autoscaling = true;
funValCheck = true;

function [stop, info] = outfun(p, optimValues, state)
  stop = false;
  info = "";
  global shouldPlot;
  global ynglobal;
  x = optimValues.model_x;
  y = optimValues.model_y;
  observations = optimValues.observations;
  
  if shouldPlot == true
    plot99 = figure(99);
    plot(x, ynglobal, ".k;Data points;", "markersize", 10, x, y, strcat({"r;Fit at "},  asctime (localtime (time)), ";"), "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("Current refinement step");
  endif
endfunction

function saveFiles(name, x, y, observations, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id) 
  a3min = a3-da3;
  
  La    = (nu+1)/alpha;
  lm    = nu/alpha;
  kapa  = 1/nu;

  Lc    = a3*(mu+1)/beta;
  kapc  = 1/mu;
  Nm    = mu/beta;
  N     = (mu+1)/beta;

  clear output;
  
  output(1, 1) = cellstr("Parameter");
  output(1, 2) = "Value";

  output(2, 1) = cellstr("cno");
  output(2, 2) = cno;

  output(3, 1) = cellstr("mu");
  output(3, 2) = mu;

  output(4, 1) = cellstr("beta");
  output(4, 2) = beta;

  output(5, 1) = cellstr("a3");
  output(5, 2) = a3;

  output(6, 1) = cellstr("da3");
  output(6, 2) = da3;

  output(7, 1) = cellstr("a3min");
  output(7, 2) = a3min;

  output(8, 1) = cellstr("sig3");
  output(8, 2) = sig3;

  output(9, 1) = cellstr("u3");
  output(9, 2) = u3;

  output(10, 1) = cellstr("eta");
  output(10, 2) = eta;

  output(11, 1) = cellstr("nu");
  output(11, 2) = nu;

  output(12, 1) = cellstr("alpha");
  output(12, 2) = alpha;

  output(13, 1) = cellstr("sig1");
  output(13, 2) = sig1;

  output(14, 1) = cellstr("lcc");
  output(14, 2) = lcc;

  output(15, 1) = cellstr("q");
  output(15, 2) = q;

  output(16, 1) = cellstr("dan");
  output(16, 2) = dan;

  output(17, 1) = cellstr("k");
  output(17, 2) = k;

  output(18, 1) = cellstr("const1");
  output(18, 2) = const1;

  output(19, 1) = cellstr("const2");
  output(19, 2) = const2;

  output(20, 1) = cellstr("g");
  output(20, 2) = g;

  output(21, 1) = cellstr("La");
  output(21, 2) = La;

  output(22, 1) = cellstr("lm");
  output(22, 2) = lm;

  output(23, 1) = cellstr("kapa");
  output(23, 2) = kapa;
  
  output(24, 1) = cellstr("Nm");
  output(24, 2) = Nm;

  output(25, 1) = cellstr("N");
  output(25, 2) = N;

  output(26, 1) = cellstr("Lc");
  output(26, 2) = Lc;

  output(27, 1) = cellstr("kapc");
  output(27, 2) = kapc;
  
  output(1, 4) = cellstr("s");
  output(1, 5) = cellstr("iObs(s)");
  output(1, 6) = cellstr("Fit");
  
  for i=1:length(x)
    output(i+1, 4) = x(i);
    output(i+1, 5) = observations(i);
    output(i+1, 6) = y(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", name, ".csv"), output, ";");

  filename = strcat(fitPath, "/output_", id, "_", name, ".txt");
  fid = fopen (filename, "w");
  fputs (fid, strcat("cno    = ",  num2str(cno), "\nmu     = ", num2str(mu), "\nbeta   = ", num2str(beta), "\na3     = ", num2str(a3), "\nda3    = ", num2str(da3), "\nsig3   = ", num2str(sig3), "\neta    = ", num2str(eta), "\nnu     = ", num2str(nu), "\nalpha  = ", num2str(alpha), "\nsig1   = ", num2str(sig1), "\nlcc    = ", num2str(lcc), "\nq      = ", num2str(q), "\ndan    = ", num2str(dan), "\nk      = ", num2str(k), "\nconst1 = ", num2str(const1), "\nconst2 = ", num2str(const2), "\ng      = ", num2str(g)));
  fclose (fid);
endfunction

#Additional parameters for fit tweaking
options1.AutoScaling = autoscaling;
options1.FunValCheck = funValCheck;
options1.MaxIter = maxIter;
options1.TolFun = tolFun;
options1.lbound = lb1;
options1.ubound = ub1;
options1.weights = wtWeight;
options1.TypicalX = typicalX;
options1.user_interaction = @outfun;

#options1.inequc = ; Additional constraints: Further inequality constraints. Cell-array containing up to four entries, two entries for linear inequality constraints and/or one or two entries for general inequality constraints. Either linear or general constraints may be the first entries, but the two entries for linear constraints must be adjacent and, if two entries are given for general constraints, they also must be adjacent. The two entries for linear constraints are a matrix (say m) and a vector (say v), specifying linear inequality constraints of the form m.' * parameters + v >= 0. The first entry for general constraints must be a differentiable column-vector valued function (say h), specifying general inequality constraints of the form h (p[, idx]) >= 0; p is the column vector of optimized parameters and the optional argument idx is a logical index. h has to return the values of all constraints if idx is not given. It may choose to return only the indexed constraints if idx is given (so computation of the other constraints can be spared); in this case, the additional setting f_inequc_idx has to be set to true. In gradient determination, this function may be called with an informational third argument, whose content depends on the function for gradient determination. If a second entry for general inequality constraints is given, it must be a function computing the jacobian of the constraints with respect to the parameters. For this function, the description of the setting dfdp, see dfdp, applies, with 2 exceptions: 1) it is called with 3 arguments since it has an additional argument idx, a logical index, at second position, indicating which rows of the jacobian must be returned (if the function chooses to return only indexed rows, the additional setting df_inequc_idx has to be set to true). 2) the default jacobian function calls h with 3 arguments, since the argument idx is also supplied. Note that specifying linear constraints as general constraints will generally waste performance, even if further, non-linear, general constraints are also specified. 
#options1.equc = ; Equality constraints. Specified the same way as inequality constraints (see inequc above). 
#options1.dfdp = ; Function computing the Jacobian of the residuals with respect to the parameters, assuming residuals are reshaped to a column vector. Default: real finite differences. Will be called with the column vector of parameters and an informational structure as arguments. If dfdp was specified by the user, the informational structure has the fields f: value of residuals for current parameters, reshaped to a column vector, fixed: logical vector indicating which parameters are not optimized, so these partial derivatives need not be computed and can be set to zero, diffp, diff_onesided, lbound, ubound: identical to the user settings of this name, plabels: 1-dimensional cell-array of column-cell-arrays, each column with labels for all parameters; the first column contains the numerical indices of the parameters; the second and third columns, present for structure based parameter handling, see Parameter structures, contain the names of the parameters and the subindices of the parameters, see Non-scalar parameters, respectively. The default jacobian function will call the model function with the second argument set with fields f: as the f passed to the jacobian function, plabels: cell-array of 1x1 cell-arrays with the entries of the column-cell-arrays of plabels as passed to the jacobian function corresponding to current parameter, side: 0 for one-sided interval, 1 or 2, respectively, for the sides of a two-sided interval, and parallel: logical scalar indicating parallel computation of partial derivatives. This information can be useful if the model function can omit some computations depending on the currently computed partial derivative. 

options2.AutoScaling = autoscaling;
options2.FunValCheck = funValCheck;
options2.MaxIter = maxIter;
options2.TolFun = tolFun;
options2.lbound = lb2;
options2.ubound = ub2;
options2.weights = wtWeight;
options2.TypicalX = typicalX;
options2.user_interaction = @outfun;

options3.AutoScaling = autoscaling;
options3.FunValCheck = funValCheck;
options3.MaxIter = maxIter;
options3.TolFun = tolFun;
options3.lbound = lb3;
options3.ubound = ub3;
options3.weights = wtWeight;
options3.TypicalX = typicalX;
options3.user_interaction = @outfun;

options4.AutoScaling = autoscaling;
options4.FunValCheck = funValCheck;
options4.MaxIter = maxIter;
options4.TolFun = tolFun;
options4.lbound = lb4;
options4.ubound = ub4;
options4.weights = wtWeight;
options4.TypicalX = typicalX;
options4.user_interaction = @outfun;

options5.AutoScaling = autoscaling;
options5.FunValCheck = funValCheck;
options5.MaxIter = maxIter;
options5.TolFun = tolFun;
options5.lbound = lb5;
options5.ubound = ub5;
options5.weights = wtNormal;
options5.TypicalX = typicalX;
options5.user_interaction = @outfun;

#Load optim package
pkg load optim;

#Statistics
settings.ret_dfdp = false;
settings.ret_covd = false;
settings.ret_covp = true;
settings.ret_corp = true;
settings.objf_type = "wls";
#settings.dfdp = dfdp;

#Plot of start values
yStart = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
if shouldPlot == true
	plot0 = figure(100);
	plot(x, yn, ".k;Data points;", "markersize", 10, x, yStart, "r;Startdata;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("0 - Start");
endif

copyfile(filename, strcat(fitPath, '/', filename));

if plotOnly == true
  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yStart(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", "0-plotOnly", ".csv"), output, ";");

  saveFiles("0-plotOnly", x, yStart, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);
  
  saveas(plot0, (strcat(fitPath, "/0_plotOnly_", id, ".png")));
else
  #Text output 
  function y = flag(flag)
    if flag == -1
      y = "Canceled.";
    endif
    if flag == 0
      y = "Maximum id of iterations.";
    endif
    if flag == 2
      y = "Change in refinement parameters too small.";
    endif
    if flag == 3
      y = "Change in calculated function too small.";
    endif
  endfunction

  #Normalization
  fun3 = @(a, x) (fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, a(1), cH, cN, cO, cS, a(2), a(3), a(4), a(5), useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, a(6), useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc));
  "\n\n\n Normalization"
  settings.weights = options3.weights;
  function result3 = result3fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param3, yn, settings, errorCount = 1)
   maxerrorCount = 10;
    result3.covp = 100*eye(length(param3));
    try
      result3 = curvefit_stat(@(a) (fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, a(1), cH, cN, cO, cS, a(2), a(3), a(4), a(5), useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, a(6), useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc)), param3, x, yn, settings);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        result3 = result3fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param3, yn, settings, errorCount + 1);
      endif
    end_try_catch
  endfunction
    
  function [param3, f3, cvg3, outp3, result3] = fit3(fun3, paramn3, x, yn, options3, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount = 1)
    maxerrorCount = 10;
    try
      [param3, f3, cvg3, outp3] = nonlin_curvefit(fun3, paramn3, x, yn, options3);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        [param3, f3, cvg3, outp3, result3] = fit3(fun3, paramn3, x, yn, options3, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount + 1);
      else
        error("Too many errors when refine")
      endif
    end_try_catch
    result3 = result3fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param3, yn, settings);
  endfunction

  [param3, f3, cvg3, outp3, result3] = fit3(fun3, paramn3, x, yn, options3, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);
  paramn3 = param3;

  q      = param3(1)
  dan    = param3(2)
  k      = param3(3)
  const1 = param3(4)
  const2 = param3(5)
  g      = param3(6)

  convergence3 = flag(cvg3)
  outp3 = outp3;
  for i=1:length(param3)
    stdabw3(i) = sqrt(result3.covp(i, i));
    mat3(i, 1) = param3(i);
    mat3(i, 2) = 1*stdabw3(i);
    mat3(i, 3) = 2*stdabw3(i);
    mat3(i, 4) = 3*stdabw3(i);
  endfor
  mat3	= mat3
  stdabw3	= stdabw3';

  yFit3 = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
 
  if shouldPlot == true
    plot3 = figure(3);
    plot(x, yn, ".k;Data points;", "markersize", 10, x, yFit3, "r;Fit3;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("3 - Normalization");
  endif

  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yFit3(i);
  endfor
    
  cell2csv(strcat(fitPath, "/output_", id, "_", "3-normalization", ".csv"), output, ";");

  saveFiles("3-normalization", x, yFit3, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);

  #Interlayer
  fun1 = @(a, x) (fun(cno, a(1), a(2), a(3), a(4), a(5), u3, a(6), nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc));
  "\n\n\n Interlayer"
  settings.weights = options1.weights;
  function result1 = result1fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param1, yn, settings, errorCount = 1)
   maxerrorCount = 10;
    result1.covp = 100*eye(length(param1));
    try
      result1 = curvefit_stat(@(a) (fun(cno, a(1), a(2), a(3), a(4), a(5), u3, a(6), nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc)), param1, x, yn, settings);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        result1 = result1fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param1, yn, settings, errorCount + 1);
      endif
    end_try_catch
  endfunction
    
  function [param1, f1, cvg1, outp1, result1] = fit1(fun1, paramn1, x, yn, options1, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount = 1)
    maxerrorCount = 10;
    try
      [param1, f1, cvg1, outp1] = nonlin_curvefit(fun1, paramn1, x, yn, options1);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        [param1, f1, cvg1, outp1, result1] = fit1(fun1, paramn1, x, yn, options1, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount + 1);
      else
        error("Too many errors when refine")
      endif
    end_try_catch
    result1 = result1fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param1, yn, settings);
  endfunction

  [param1, f1, cvg1, outp1, result1] = fit1(fun1, paramn1, x, yn, options1, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);
  paramn1 = param1;

  mu   = param1(1)
  beta = param1(2)
  a3   = param1(3)
  da3  = param1(4)
  sig3 = param1(5)
  eta  = param1(6)

  convergence1 = flag(cvg1)
  outp1 = outp1;
  for i=1:length(param1)
    stdabw1(i) = sqrt(result1.covp(i, i));
    mat1(i, 1) = param1(i);
    mat1(i, 2) = 1*stdabw1(i);
    mat1(i, 3) = 2*stdabw1(i);
    mat1(i, 4) = 3*stdabw1(i);
  endfor
  stdabw1 = stdabw1';
  mat1 = mat1

  yFit1 = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
  
  if shouldPlot == true
    plot1 = figure(1);
    plot(x, yn, ".k;Data points;", "markersize", 10, x, yFit1, "r;Fit1;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("1 - Interlayer");
  endif

  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yFit1(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", "1-interlayer", ".csv"), output, ";");

  saveFiles("1-interlayer", x, yFit1, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);

  #Intralayer
  fun2 = @(a, x) (fun(cno, mu, beta, a3, da3, sig3, u3, eta, a(1), a(2), lcc, a(3), q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc));
  "\n\n\n Intralayer"
  settings.weights = options2.weights;
  function result2 = result2fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param2, yn, settings, errorCount = 1)
   maxerrorCount = 10;
    result2.covp = 100*eye(length(param2));
    try
      result2 = curvefit_stat(@(a) (fun(cno, mu, beta, a3, da3, sig3, u3, eta, a(1), a(2), lcc, a(3), q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc)), param2, x, yn, settings);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        result2 = result2fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param2, yn, settings, errorCount + 1);
      endif
    end_try_catch
  endfunction
    
  function [param2, f2, cvg2, outp2, result2] = fit2(fun2, paramn2, x, yn, options2, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount = 1)
    maxerrorCount = 10;
    try
      [param2, f2, cvg2, outp2] = nonlin_curvefit(fun2, paramn2, x, yn, options2);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        [param2, f2, cvg2, outp2, result2] = fit2(fun2, paramn2, x, yn, options2, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount + 1);
      else
        error("Too many errors when refine")
      endif
    end_try_catch
    result2 = result2fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param2, yn, settings);
   endfunction
   
  [param2, f2, cvg2, outp2, result2] = fit2(fun2, paramn2, x, yn, options2, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);
  paramn2 = param2;

  nu    = param2(1)
  alpha = param2(2)
  sig1  = param2(3)

  convergence2 = flag(cvg2)
  outp2 = outp2;
  for i=1:length(param2)
    stdabw2(i) = sqrt(result2.covp(i, i));
    mat2(i, 1) = param2(i);
    mat2(i, 2) = 1*stdabw2(i);
    mat2(i, 3) = 2*stdabw2(i);
    mat2(i, 4) = 3*stdabw2(i);
  endfor
  stdabw2 = stdabw2';
  mat2 = mat2

  yFit2 = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
  
  if shouldPlot == true
    plot2 = figure(2);
    plot(x, yn, ".k;Data points;", "markersize", 10, x, yFit2, "r;Fit2;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("2 - Intralayer");
  endif

  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yFit2(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", "2-intralayer", ".csv"), output, ";");

  saveFiles("2-intralayer", x, yFit2, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);

  #Normalisierung
  "\n\n\n Normalization"
  [param3, f3, cvg3, outp3, result3] = fit3(fun3, paramn3, x, yn, options3, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);
  paramn3 = param3;

  q      = param3(1)
  dan    = param3(2)
  k      = param3(3)
  const1 = param3(4)
  const2 = param3(5)
  g      = param3(6)

  convergence3 = flag(cvg3)
  outp3 = outp3;
  for i=1:length(param3)
    stdabw3(i) = sqrt(result3.covp(i, i));
    mat3(i, 1) = param3(i);
    mat3(i, 2) = 1*stdabw3(i);
    mat3(i, 3) = 2*stdabw3(i);
    mat3(i, 4) = 3*stdabw3(i);
  endfor
  mat3 = mat3
  stdabw3 = stdabw3';

  yFit3 = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
  
  if shouldPlot == true
    plot3 = figure(3);
    plot(x, yn, ".k;Data points;", "markersize", 10, x, yFit3, "r;Fit3;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("3 - Normalization");
  endif

  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yFit3(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", "3-normalization", ".csv"), output, ";");

  saveFiles("3-normalization", x, yFit3, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);

  #lcc
  fun4 = @(a, x) (fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, a(1), sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc));
  "\n\n\n lcc"
  settings.weights = options4.weights;
  function result4 = result4fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param4, yn, settings, errorCount = 1)
   maxerrorCount = 10;
    result4.covp = 100*eye(length(param4));
    try
      result4 = curvefit_stat(@(a) (fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, a(1), sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc)), param4, x, yn, settings);     
    catch
      lasterror.message
      if errorCount < maxerrorCount
        result4 = result4fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param4, yn, settings, errorCount + 1);
      endif
    end_try_catch
  endfunction
    
  function [param4, f4, cvg4, outp4, result4] = fit4(fun4, paramn4, x, yn, options4, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount = 1)
    maxerrorCount = 10;
    try
      [param4, f4, cvg4, outp4] = nonlin_curvefit(fun4, paramn4, x, yn, options4);
     catch
      lasterror.message
      if errorCount < maxerrorCount
        [param4, f4, cvg4, outp4, result4] = fit4(fun4, paramn4, x, yn, options4, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, errorCount + 1);
      else
        error("Too many errors when refine")
      endif
    end_try_catch
    result4 = result4fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, x, param4, yn, settings);
  endfunction

  [param4, f4, cvg4, outp4, result4] = fit4(fun4, paramn4, x, yn, options4, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);
  paramn4 = param4;

  lcc = param4(1)

  convergence4 = flag(cvg4)
  outp4 = outp4;
  for i=1:length(param4)
    stdabw4(i) = sqrt(result4.covp(i, i));
    mat4(i, 1) = param4(i);
    mat4(i, 2) = 1*stdabw4(i);
    mat4(i, 3) = 2*stdabw4(i);
    mat4(i, 4) = 3*stdabw4(i);
  endfor
  stdabw4 = stdabw4';
  mat4 = mat4

  yFit4 = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
  
  if shouldPlot == true
    plot4 = figure(4);
    plot(x, yn, ".k;Data points;", "markersize", 10, x, yFit4, "r;Fit4;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("4 - lcc");
  endif

  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yFit4(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", "4-lcc", ".csv"), output, ";");

  saveFiles("4-lcc", x, yFit4, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);

  #All
  fun5 = @(a, x) (fun(cno, a(1), a(2), a(3), a(4), a(5), u3, a(6), a(7), a(8), a(9), a(10), a(11), cH, cN, cO, cS, a(12), a(13), a(14), a(15), useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, a(16), useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc));
  "\n\n\n All"
  settings.weights = options5.weights;
  paramn5 = [mu; beta; a3; da3; sig3; eta; nu; alpha; lcc; sig1; q; dan; k; const1; const2; g];
  function result5 = result5fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, coh, inc, x, param5, yn, settings, errorCount = 1)
   maxerrorCount = 10;
    result5.covp = 100*eye(length(param5));
    try
      result5 = curvefit_stat(@(a) (fun(cno, a(1), a(2), a(3), a(4), a(5), u3, a(6), a(7), a(8), a(9), a(10), a(11), cH, cN, cO, cS, a(12), a(13), a(14), a(15), useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, a(16), useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc)), param5, x, yn, settings);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        result5 = result5fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, coh, inc, x, param5, yn, settings, errorCount + 1);
      endif
    end_try_catch
  endfunction
    
  function [param5, f5, cvg5, outp5, result5] = fit5(fun5, paramn5, x, yn, options5, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, coh, inc, errorCount = 1)
    maxerrorCount = 10;
    try
      [param5, f5, cvg5, outp5] = nonlin_curvefit(fun5, paramn5, x, yn, options5);
    catch
      lasterror.message
      if errorCount < maxerrorCount
        [param5, f5, cvg5, outp5, result5] = fit5(fun5, paramn5, x, yn, options5, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, coh, inc, errorCount + 1);
      else
        error("Too many errors when refine")
      endif
    end_try_catch
    result5 = result5fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, coh, inc, x, param5, yn, settings);
  endfunction

  [param5, f5, cvg5, outp5, result5] = fit5(fun5, paramn5, x, yn, options5, settings, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, coh, inc);
  paramn5 = param5;

  mu     = param5(1)
  beta   = param5(2)
  a3     = param5(3)
  da3    = param5(4)
  sig3   = param5(5)
  eta    = param5(6)
  nu     = param5(7)
  alpha  = param5(8)
  lcc    = param5(9)
  sig1   = param5(10)
  q      = param5(11)
  dan    = param5(12)
  k      = param5(13)
  const1 = param5(14)
  const2 = param5(15)
  g      = param5(16)

  convergence5 = flag(cvg5)
  outp5 = outp5;
  for i=1:length(param5)
    stdabw5(i) = sqrt(result5.covp(i, i));
    mat5(i, 1) = param5(i);
    mat5(i, 2) = 1*stdabw5(i);
    mat5(i, 3) = 2*stdabw5(i);
    mat5(i, 4) = 3*stdabw5(i);
  endfor
  stdabw5 = stdabw5';
  mat5 = mat5

  yFit5 = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
   
  if shouldPlot == true
    plot5 = figure(5);
    plot(x, yn, ".k;Data points;", "markersize", 10, x, yFit5, "r;Fit5;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("5 - All");
  endif

  output(1, 1) = cellstr("s");
  output(1, 2) = cellstr("iObs(s)");
  for i=1:length(x)
    output(i+1, 1) = x(i);
    output(i+1, 2) = yFit5(i);
  endfor
  
  cell2csv(strcat(fitPath, "/output_", id, "_", "5-all", ".csv"), output, ";");

  saveFiles("5-all", x, yFit5, yn, cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, s, coh, inc, fitPath, id);

  paramFinish = [cno; mu; beta; a3; da3; sig3; u3; eta; nu; alpha; sig1; lcc; q; dan; k; const1; const2; g];

  yFit = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
  yFitPlot = fun(cno, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, lcc, sig1, q, cH, cN, cO, cS, dan, k, const1, const2, useQ, b, useA, density, sampleThickness, transmission, absorptionCorrection, useP, polarizedBeam, polarizationDegree, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, radiation, wavelength, x, coh, inc);
  dy = ((yn-yFit)./yn)';
  dySumme = sum(dy)

  x0 = [min(x); max(x)];
  y0 = [0; 0];
  yFitLog = log10(yFit);
  yFitLogPlot = log10(yFitPlot);

  ynLog = log10(yn);
  dyLog = (ynLog-yFitLog)./ynLog;
  dyLogSumme = sum(dyLog)

  x0 = [min(x); max(x)];
  y0 = [0; 0];

  if shouldPlot == true
    plot7 = figure(8);
    plot(x0, y0, "k;Zero line;", "LineWidth", 3, x, dy, ".r;errorCount;");
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("7 - Error");
    plot8 = figure(8);
    plot(x, ynLog, ".k;Data points;", "markersize", 10, x, yFitLogPlot, "r;Fit;", "LineWidth", 3);
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("8 - Log");
    plot9 = figure(9);
    plot(x0, y0, "k;Zero line;", "LineWidth", 3, x, dyLog, ".r;errorCount;");
	xlabel ("Scattering vector s / A^-^1");
	ylabel ("Intensity I");
	title ("9 - Log error");
  endif

  rQuadratFit = regress(yFit', yn');

  chiQuadratFit = 0;
  for i=1:length(yn)
    chiQuadratFit = chiQuadratFit + log10(wtNormal(i)*(yn(i)-yFit(i))^2);
  endfor

  chiQuadratFit = chiQuadratFit

  if (rQuadratFit < 1)
    rQuadratFit = rQuadratFit
  else
    rQuadratFit = 1/rQuadratFit
endif

  if shouldPlot == true
    saveas(plot0, (strcat(fitPath, "/0_start_", id, ".png")));
    saveas(plot1, (strcat(fitPath, "/1_interlayer_", id, ".png")));
    saveas(plot2, (strcat(fitPath, "/2_intralayer_", id, ".png")));
    saveas(plot3, (strcat(fitPath, "/3_normalisierung_", id, ".png")));
    saveas(plot4, (strcat(fitPath, "/4_lcc_", id, ".png")));
    saveas(plot5, (strcat(fitPath, "/5_all_", id, ".png")));
    saveas(plot7, (strcat(fitPath, "/7_errorCount_", id, ".png")));
    saveas(plot8, (strcat(fitPath, "/8_allLog_", id, ".png")));
    saveas(plot9, (strcat(fitPath, "/9_errorCountLog_", id, ".png")));
  endif

  save([fitPath, strcat("/data_", id, ".mat")],
  "tolFun", "maxIter", "paramn", "paramFinish",
  "param1", "convergence1", "outp1", "result1", "stdabw1", "mat1",
  "param2", "convergence2", "outp2", "result2", "stdabw2", "mat2",
  "param3", "convergence3", "outp3", "result3", "stdabw3", "mat3",
  "param4", "convergence4", "outp4", "result4", "stdabw4", "mat4",
  "param5", "convergence5", "outp5", "result5", "stdabw5", "mat5",
  "rQuadratFit", "chiQuadratFit");

  "\n\n\n Refined values"
  cno      = cno
  cno1     = 0;
  cno2     = 2*cno1;
  cno3     = 3*cno1;
  mu       = mu
  mu1      = mat5(1, 2);
  mu2      = 2*mu1;
  mu3      = 3*mu1;
  beta     = beta
  beta1    = mat5(2, 2);
  beta2    = 2*beta1;
  beta3    = 3*beta1;
  a3       = a3
  a31      = mat5(3, 2);
  a32      = 2*a31;
  a33      = 3*a31;
  da3      = da3
  da31     = mat5(4, 2);
  da32     = 2*da31;
  da33     = 3*da31;
  a3min    = a3-da3
  a3min1   = a31+da31;
  a3min2   = 2*a3min1;
  a3min3   = 3*a3min1;
  sig3     = sig3
  sig31    = mat5(5, 2);
  sig32    = 2*sig31;
  sig33    = 3*sig31;
  u3       = u3;
  u31      = 0;
  u32      = 2*u31;
  u33      = 3*u31;
  eta      = eta
  eta1     = mat5(6, 2);
  eta2     = 2*eta1;
  eta3     = 3*eta1;
  nu       = nu
  nu1      = 0;
  nu2      = 2*nu1;
  nu3      = 3*nu1;
  alpha    = alpha
  alpha1   = mat5(7, 2);
  alpha2   = 2*alpha1;
  alpha3   = 3*alpha1;
  sig1     = sig1
  sig11    = mat5(8, 2);
  sig12    = 2*sig11;
  sig13    = 3*sig11;
  lcc      = lcc
  lcc1     = mat5(9, 2);
  lcc2     = 2*lcc1;
  lcc3     = 3*lcc1;
  q        = q
  q1       = mat5(10, 2);
  q2       = 2*q1;
  q3       = 3*q1;
  dan      = dan
  dan1     = mat5(11, 2);
  dan2     = 2*dan1;
  dan3     = 3*dan1;
  k        = k
  k1       = mat5(12, 2);
  k2       = 2*k1;
  k3       = 3*k1;
  const1   = const1
  const11  = mat5(13, 2);
  const12  = 2*const11;
  const13  = 3*const11;
  const2   = const2
  const21  = mat5(14, 2);
  const22  = 2*const21;
  const23  = 3*const21;
  g        = g
  g1       = mat5(15, 2);
  g2       = 2*g1;
  g3       = 3*g1;

  "\n\n\n Intralayer"
  La       = (nu+1)/alpha
  La1      = (nu1/alpha+(1+nu)/alpha^2*alpha1)/2;
  La2      = 2*La1;
  La3      = 3*La1;
  lm       = nu/alpha
  lm1      = nu1/alpha+nu/alpha^2*alpha1;
  lm2      = 2*lm1;
  lm3      = 3*lm1;
  kapa     = 1/nu
  kapa1    = nu1/nu^2;
  kapa2    = 2*kapa1;
  kapa3    = 3*kapa1;
  lcc      = lcc
  lcc1     = lcc1;
  lcc2     = 2*lcc1;
  lcc3     = 3*lcc1;
  sig1     = sig1
  sig11    = sig11;
  sig12    = 2*sig11;
  sig13    = 3*sig11;

  "\n\n\n Interlayer"
  Lc       = a3*(mu+1)/beta
  Lc1      = (mu+1)/beta*a31+(a3*beta1/beta+a3*(1+mu)/beta^2*beta1)/2;
  Lc2      = 2*Lc1;
  Lc3      = 3*Lc1;
  kapc     = 1/mu
  kapc1    = mu1/mu^2;
  kapc2    = 2*kapc1;
  kapc3    = 3*kapc1;
  a3       = a3
  a31      = a31;
  a32      = 2*a31;
  a33      = 3*a31;
  a3min    = a3min
  a3min1   = a3min1;
  a3min2   = 2*a3min1;
  a3min3   = 3*a3min1;
  sig3	 = sig3
  sig31    = sig31;
  sig32    = 2*sig31;
  sig33    = 3*sig31;
  Nm       = mu/beta
  Nm1      = mu1/beta+mu/beta^2*beta1;
  Nm2      = 2*Nm1;
  Nm3      = 3*Nm1;
  N        = (mu+1)/beta
  N1       = (mu1/beta+(1+mu)/beta^2*beta1)/2;
  N2       = 2*N1;
  N3       = 3*N1;
  u3       = u3
  u31      = 0;
  u32      = 2*u31;
  u33      = 3*u31;
  eta      = eta
  eta1     = eta1;
  eta2     = 2*eta1;
  eta3     = 3*eta1;
  q        = q
  q1       = q1;
  q2       = 2*q1;
  q3       = 3*q1;
  dan      = dan
  dan1     = 0;
  dan2     = 2*dan1;
  dan3     = 3*dan1;

  output(1, 1) = cellstr("Parameter");
  output(1, 2) = "Value";
  output(1, 3) = "1-sigma errorCount";
  output(1, 4) = "2-sigma errorCount";
  output(1, 5) = "3-sigma errorCount";

  output(2, 1) = cellstr("cno");
  output(2, 2) = cno;
  output(2, 3) = cno1;
  output(2, 4) = cno2;
  output(2, 5) = cno3;

  output(3, 1) = cellstr("mu");
  output(3, 2) = mu;
  output(3, 3) = mu1;
  output(3, 4) = mu2;
  output(3, 5) = mu3;

  output(4, 1) = cellstr("beta");
  output(4, 2) = beta;
  output(4, 3) = beta1;
  output(4, 4) = beta2;
  output(4, 5) = beta3;

  output(5, 1) = cellstr("a3");
  output(5, 2) = a3;
  output(5, 3) = a31;
  output(5, 4) = a32;
  output(5, 5) = a33;

  output(6, 1) = cellstr("da3");
  output(6, 2) = da3;
  output(6, 3) = da31;
  output(6, 4) = da32;
  output(6, 5) = da33;

  output(7, 1) = cellstr("a3min");
  output(7, 2) = a3min;
  output(7, 3) = a3min1;
  output(7, 4) = a3min2;
  output(7, 5) = a3min3;

  output(8, 1) = cellstr("sig3");
  output(8, 2) = sig3;
  output(8, 3) = sig31;
  output(8, 4) = sig32;
  output(8, 5) = sig33;

  output(9, 1) = cellstr("u3");
  output(9, 2) = u3;
  output(9, 3) = u31;
  output(9, 4) = u32;
  output(9, 5) = u33;

  output(10, 1) = cellstr("eta");
  output(10, 2) = eta;
  output(10, 3) = eta1;
  output(10, 4) = eta2;
  output(10, 5) = eta3;

  output(11, 1) = cellstr("nu");
  output(11, 2) = nu;
  output(11, 3) = nu1;
  output(11, 4) = nu2;
  output(11, 5) = nu3;

  output(12, 1) = cellstr("alpha");
  output(12, 2) = alpha;
  output(12, 3) = alpha1;
  output(12, 4) = alpha2;
  output(12, 5) = alpha3;

  output(13, 1) = cellstr("sig1");
  output(13, 2) = sig1;
  output(13, 3) = sig11;
  output(13, 4) = sig12;
  output(13, 5) = sig13;

  output(14, 1) = cellstr("lcc");
  output(14, 2) = lcc;
  output(14, 3) = lcc1;
  output(14, 4) = lcc2;
  output(14, 5) = lcc3;

  output(15, 1) = cellstr("q");
  output(15, 2) = q;
  output(15, 3) = q1;
  output(15, 4) = q2;
  output(15, 5) = q3;

  output(16, 1) = cellstr("dan");
  output(16, 2) = dan;
  output(16, 3) = dan1;
  output(16, 4) = dan2;
  output(16, 5) = dan3;

  output(17, 1) = cellstr("k");
  output(17, 2) = k;
  output(17, 3) = k1;
  output(17, 4) = k2;
  output(17, 5) = k3;

  output(18, 1) = cellstr("const1");
  output(18, 2) = const1;
  output(18, 3) = const11;
  output(18, 4) = const12;
  output(18, 5) = const13;

  output(19, 1) = cellstr("const2");
  output(19, 2) = const2;
  output(19, 3) = const21;
  output(19, 4) = const22;
  output(19, 5) = const23;

  output(20, 1) = cellstr("g");
  output(20, 2) = g;
  output(20, 3) = g1;
  output(20, 4) = g2;
  output(20, 5) = g3;

  output(21, 1) = cellstr("La");
  output(21, 2) = La;
  output(21, 3) = La1;
  output(21, 4) = La2;
  output(21, 5) = La3;

  output(22, 1) = cellstr("lm");
  output(22, 2) = lm;
  output(22, 3) = lm1;
  output(22, 4) = lm2;
  output(22, 5) = lm3;

  output(23, 1) = cellstr("kapa");
  output(23, 2) = kapa;
  output(23, 3) = kapa1;
  output(23, 4) = kapa2;
  output(23, 5) = kapa3;

  output(24, 1) = cellstr("Nm");
  output(24, 2) = Nm;
  output(24, 3) = Nm1;
  output(24, 4) = Nm2;
  output(24, 5) = Nm3;

  output(25, 1) = cellstr("N");
  output(25, 2) = N;
  output(25, 3) = N1;
  output(25, 4) = N2;
  output(25, 5) = N3;

  output(26, 1) = cellstr("Lc");
  output(26, 2) = Lc;
  output(26, 3) = Lc1;
  output(26, 4) = Lc2;
  output(26, 5) = Lc3;

  output(27, 1) = cellstr("kapc");
  output(27, 2) = kapc;
  output(27, 3) = kapc1;
  output(27, 4) = kapc2;
  output(27, 5) = kapc3;

  output(1, 7) = "s";
  output(1, 8) = "I";
  output(1, 9) = "Ifit";
  output(1, 10) = "abwAbs";
  output(1, 11) = "abwRel";

  abwAbs = (yn-yFit);
  abwRel = dy;
  for i=1:length(x)
    output(i+1, 7) = x(i);
    output(i+1, 8) = yn(i);
    output(i+1, 9) = yFitPlot(i);
    output(i+1, 10) = abwAbs(i);
    output(i+1, 11) = abwRel(i);
  endfor

  output = output;

  cell2csv(strcat(fitPath, "/output_", id, ".csv"), output, ";");

  filename = strcat(fitPath, "/output_", id, ".txt");
  fid = fopen (filename, "w");
  fputs (fid, strcat("cno		= ",	num2str(cno), "\nmu		 = ", num2str(mu), "\nbeta	 = ", num2str(beta), "\na3		 = ", num2str(a3), "\nda3		= ", num2str(da3), "\nsig3	 = ", num2str(sig3), "\neta		= ", num2str(eta), "\nnu		 = ", num2str(nu), "\nalpha	= ", num2str(alpha), "\nsig1	 = ", num2str(sig1), "\nlcc		= ", num2str(lcc), "\nq			= ", num2str(q), "\ndan		= ", num2str(dan), "\nk			= ", num2str(k), "\nconst1 = ", num2str(const1), "\nconst2 = ", num2str(const2), "\ng			= ", num2str(g)));
  fclose (fid);
endif

timeEnd = time();

duration = timeEnd-timeStart

#Copyright (C) 2021 Oliver Osswald