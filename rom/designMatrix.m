function [Phi] = designMatrix(phi, fineCond, domainf, domainc)
%Construct design matrix Phi from basis functions phi and input conductivity data x.

%vector E gives the coarse element a fine element belongs to
[E] = get_coarse_el([domainf.nElX, domainf.nElY], [domainc.nElX, domainc.nElY], 1:domainf.nEl);

%inputs belonging to same coarse element are in the same column of xk. They are ordered in
%x-direction.
lambdak = zeros(domainf.nEl/domainc.nEl, domainc.nEl);
for i = 1:domainc.nEl
    lambdak(:, i) = fineCond(E == i);
end

%construct design matrix Phi
Phi = zeros(domainc.nEl, size(phi, 1));
for i = 1:domainc.nEl
    for j = 1:size(phi, 1)
        Phi(i, j) = phi{j}(lambdak(:, i));
    end
end

end

