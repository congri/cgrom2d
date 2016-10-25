function [PhiArray] = designMatrix(phi, domainf, domainc)
%Construct design matrix Phi from basis functions phi and input conductivity data x.

%load fineCond from disc
load('./data/fineData/fineData', 'cond');

%vector E gives the coarse element a fine element belongs to
[E] = get_coarse_el([domainf.nElX, domainf.nElY], [domainc.nElX, domainc.nElY], 1:domainf.nEl);

PhiArray = zeros(domainc.nEl, size(phi, 1), size(cond, 1));
for s = 1:size(cond, 1)
    %inputs belonging to same coarse element are in the same column of xk. They are ordered in
    %x-direction.
    lambdak = zeros(domainf.nEl/domainc.nEl, domainc.nEl);
    for i = 1:domainc.nEl
        lambdak(:, i) = cond(s, E == i);
    end
    
    %construct design matrix Phi
    for i = 1:domainc.nEl
        for j = 1:size(phi, 1)
            PhiArray(i, j, s) = phi{j}(lambdak(:, i));
        end
    end
    
end

end

