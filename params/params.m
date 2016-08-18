%main parameter file for 2d coarse-graining

%initialize domain
nF = 10;
nC = 5;
%conductivity in unit cell
lambdaUnit = [1 3; 2 4];
assert(mod(nF, 2) == 0, 'error: nF not divisible by 2')
dFine = Domain(nF, nF, 1, 1);
dCoarse = Domain(nC, nC, 1, 1);
%boundary conditions
boundaryConditions;
dFine = setNodalCoordinates(dFine, physical);
dFine = setBvec(dFine);

dCoarse = setNodalCoordinates(dCoarse, physical);
dCoarse = setBvec(dCoarse);

%physical params, boundary conditions etc.
physC = physicalParams(dCoarse, physical);
physF = physicalParams(dFine, physical);

%Finescale conductivity params
fineCond.mu = 3;    %mean of log of lambda
fineCond.sigma = .5; %sigma of log of lambda
fineCond.nSamples = 2;

%Define basis functions for p_c here
phi_1 = @(x) size(x, 1)/sum(1./x);
phi_2 = @(x) mean(x);
phi = {phi_1; phi_2};