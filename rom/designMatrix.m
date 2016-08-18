function [Phi] = designMatrix(phi, x, dFine, dCoarse)
%Construct design matrix Phi from basis functions phi and input conductivity data x. Keep in mind
%that nFine, nCoarse are vectors with elements in x- and y-direction

%vector E gives the coarse element a fine element belongs to
[E] = get_coarse_el([dFine.nElX, dFine.nElY], [dCoarse.nElX, dCoarse.nElY], 1:dFine.nEl);

%inputs belonging to same coarse element are in the same column of xk
xk = zeros(dFine.nEl/dCoarse.nEl, dCoarse.nEl);
for i = 1:dCoarse.nEl
    xk(:,i) = x(E == i);
end

%construct design matrix Phi
Phi = zeros(dCoarse.nEl, size(phi, 1));
for i = 1:dCoarse.nEl
    for j = 1:size(phi, 1)
        Phi(i, j) = phi{j}(xk(:,i));
    end
end

end

