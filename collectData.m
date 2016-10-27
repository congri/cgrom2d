%Script to collect data in data arrays of EM object

if ~exist(strcat('./data/', jobname), 'dir')
    mkdir(strcat('./data/', jobname));
end

%% MCMC Step width
MCMCStepWidth = zeros(1, nTrain);
filename = strcat('./data/', jobname, '/MCMCstepWidth');
for i = 1:nTrain
    if strcmp(MCMC(i).method, 'MALA')
        MCMCStepWidth(i) = MCMC(i).MALA.stepWidth;
    elseif strcmp(MCMC(i).method, 'randomWalk')
        %only valid for isotropic proposals!
        MCMCStepWidth(i) = MCMC(i).randomWalk.proposalCov(1, 1);
    elseif strcmp(MCMC(i).method, 'nonlocal')
        %do nothing; we don't use this
    else
        error('Unknown sampling method')
    end
end
save(filename, 'MCMCStepWidth', '-ascii', '-append')

%% Optimal params
%W matrix
saveW = true;
if saveW
    filename = strcat('./data/', jobname, '/Wmat');
    [rowW, colW, valW] = find(theta_cf.W);
    WArray = [rowW, colW, valW]';
    onlyFinal = true;
    if onlyFinal
        save(filename, 'WArray', '-ascii')
    else
        save(filename, 'WArray', '-ascii', '-append')
    end
end
clear rowW colW valW WArray;
%theta
filename = strcat('./data/', jobname, '/theta');
theta = theta_c.theta';
save(filename, 'theta', '-ascii', '-append');
%sigma
filename = strcat('./data/', jobname, '/sigma');
sigma = theta_c.sigma;
save(filename, 'sigma', '-ascii', '-append');
%S
saveS = true;
if saveS
    filename = strcat('./data/', jobname, '/S');
    S = theta_cf.S';
    onlyFinal = true;
    if onlyFinal
        save(filename, 'S', '-ascii');
    else
        save(filename, 'S', '-ascii', '-append');
    end
    clear S;
end
%mu
saveMu = true;
if saveMu
    mu = theta_cf.mu';
    filename = strcat('./data/', jobname, '/mu');
    onlyFinal = true;
    if onlyFinal
        save(filename, 'mu', '-ascii')
    else
        save(filename, 'mu', '-ascii', '-append')
    end
    clear mu;
end
