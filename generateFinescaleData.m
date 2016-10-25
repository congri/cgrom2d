%Script to generate and save finescale data
tic;
clear all;
addpath('./params')
addpath('./heatFEM')
addpath('./rom')
addpath('./genConductivity')

%% Temperature field and gradient generating the boundary conditions
a = [-5 3 -2 0];
Tb = @(x) a(1) + a(2)*x(1) + a(3)*x(2) + a(4)*x(1)*x(2);
qb{1} = @(x) -(a(3) + a(4)*x);      %lower bound
qb{2} = @(y) (a(2) + a(4)*y);       %right bound
qb{3} = @(x) (a(3) + a(4)*x);       %upper bound
qb{4} = @(y) -(a(2) + a(4)*y);      %left bound

%% Generate finescale domain
nf = 128;       %Should be 2^n
disp('Generate domain...')
domainf = Domain(nf, nf, 1, 1);
domainf = setBoundaries(domainf, 2:(4*nf), Tb, qb);           %ATTENTION: natural nodes have to be set manually
                                                                %and consistently in domainc and domainf
domainf = setNodalCoordinates(domainf);
domainf = setBvec(domainf);
domainf = setHeatSource(domainf, zeros(domainf.nEl, 1));

%% Finescale data params
t = toc
disp('Setting up finescale data parameters...')
fineData.genData = true;    %generate new dataset?
fineData.dist = 'correlated_binary';   %uniform, gaussian or binary (dist of log conductivity)
fineData.nSamples = 3;
if strcmp(fineData.dist, 'gaussian')
    fineData.mu = 1.2;      %mean of log of lambda
    fineData.sigma = .3;    %sigma of log of lambda
elseif (strcmp(fineData.dist, 'uniform') || strcmp(fineData.dist, 'binary')...
        || strcmp(fineData.dist, 'predefined_binary') || strcmp(fineData.dist, 'correlated_binary'))
    %for uniform & binary
    fineData.lo = 1;    %upper and lower bounds on conductivity lambda
    fineData.up = 1e1;
    contrast = fineData.up/fineData.lo;
    %for binary
    if (strcmp(fineData.dist, 'binary') || strcmp(fineData.dist, 'correlated_binary'))
        fineData.p_lo = 0;
        if strcmp(fineData.dist, 'correlated_binary')
            fineData.lx = 6*domainf.lElX;
            fineData.ly = 6*domainf.lElY;
            fineData.sigma_f2 = 1; %has this parameter any impact?
        end
    end
else
    error('unknown fineCond distribution');
end



%% Generate finescale data
t = toc
genData(domainf, fineData);
t = toc

disp('Finescale data generated. Saving and quitting...')
% save data
save('./data/fineData/fineDataParams', 'domainf', 'fineData', 'Tb', 'qb');

time_for_data_generation = toc