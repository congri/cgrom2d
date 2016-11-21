function [PhiArray] = designMatrix(phi, domainf, domainc, Tffile, nTrain_lo, nTrain_up)
%Construct design matrix Phi from basis functions phi and input conductivity data x.

tic
disp('Compute design matrices Phi...')
%load fineCond from disc
cond = Tffile.cond(nTrain_lo:nTrain_up, :);        %Finescale conductivities - load partially to save memory

%vector E gives the coarse element a fine element belongs to
[E] = get_coarse_el([domainf.nElX, domainf.nElY], [domainc.nElX, domainc.nElY], 1:domainf.nEl);

PhiArray = zeros(domainc.nEl, numel(phi), size(cond, 1));
parfor s = 1:size(cond, 1)
    %inputs belonging to same coarse element are in the same column of xk. They are ordered in
    %x-direction.
    PhiCell{s} = zeros(domainc.nEl, numel(phi));
    lambdak = zeros(domainf.nEl/domainc.nEl, domainc.nEl);
    for i = 1:domainc.nEl
        lambdak(:, i) = cond(s, E == i);
    end
    
    %construct design matrix Phi
    for i = 1:domainc.nEl
        for j = 1:numel(phi)
            PhiCell{s}(i, j) = phi{j}(lambdak(:, i));
        end
    end  
end

for i = 1:size(cond, 1)
    PhiArray(:, :, i) = PhiCell{i};
end
disp('done')
Phi_computation_time = toc

end

