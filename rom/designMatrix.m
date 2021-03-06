function [Phi] = designMatrix(phi, x, domainf, domainc)
%Construct design matrix Phi from basis functions phi and input conductivity data x. Keep in mind
%that nFine, nCoarse are vectors with elements in x- and y-direction

%vector E gives the coarse element a fine element belongs to
[E] = get_coarse_el([domainf.nElX, domainf.nElY], [domainc.nElX, domainc.nElY], 1:domainf.nEl);

%inputs belonging to same coarse element are in the same column of xk
xk = zeros(domainf.nEl/domainc.nEl, domainc.nEl);
for i = 1:domainc.nEl
    xk(:,i) = x(E == i);
end

%construct design matrix Phi
Phi = zeros(domainc.nEl, size(phi, 1));
for i = 1:domainc.nEl
    for j = 1:size(phi, 1)
        Phi(i, j) = phi{j}(exp(xk(:, i)));
    end
end

end

