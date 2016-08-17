%main parameter file for 2d coarse-graining

%initialize domain
dFine = Domain(16, 16, 1, 1);
%boundary conditions
boundaryConditions;
dFine = setNodalCoordinates(dFine, physical);
dFine = setBvec(dFine);

%physical params, boundary conditions etc.
physicalParams;