"Start installing the optim package"
pkg install -forge optim
"optim package installed"

"Load optim package"
pkg load optim
"optim package loaded"

"Test optim package"
try
  nonlin_curvefit
catch err
  str = err(1).message;
  if regexp(str, "Invalid call to nonlin_curvefit*") >= 1
    "Test completed successfully. The package has been installed."
  else
    "Error!"
    error = err
    "If you cannot fix the error yourself, please copy the entire issue and contact Bernd Smarsly (bernd.smarsly@phys.chemie.uni-giesesn.de) or Oliver Osswald (oliver.osswald@phys.chemie. uni-giessen.de)"
  endif
end_try_catch