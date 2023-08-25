"Start installing the optim package"
pkg install -forge optim
"optim package installed"

#The parallel computing is currently not working due to several errors (may be caused by the global variables, which are not supported by the optim package. To use it, load the package and set options1.parallel_local = true; Instead of true, a specific core number can be used
#"Start installing the parallel package"
#pkg install -forge parallel
#"parallel package installed"

"Test optim package"
try
  "Load optim package"
  pkg load optim
  "optim package loaded"
catch err
  str = err(1).message;
  if str != "package optim is not installed"
    "Test completed successfully. The optim package has been installed."
  else
    "Error!"
    error = err
    "If you cannot fix the error yourself, please copy the entire issue and contact Bernd Smarsly (bernd.smarsly@phys.chemie.uni-giesesn.de) or Oliver Osswald (oliver.osswald@phys.chemie.uni-giessen.de)"
  endif
end_try_catch

#"Test parallel package"
#try
#  "Load parallel package"
#  pkg load parallel
#  "parallel package loaded"
#catch err
#  str = err(1).message;
#  if str != "package parallel is not installed"
#    "Test completed successfully. The parallel package has been installed."
#  else
#    "Error!"
#    error = err
#    "If you cannot fix the error yourself, please copy the entire issue and contact Bernd Smarsly (bernd.smarsly@phys.chemie.uni-giesesn.de) or Oliver Osswald (oliver.osswald@phys.chemie.uni-giessen.de)"
#  endif
#end_try_catch