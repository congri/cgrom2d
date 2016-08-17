%specify boundary conditions here
physical.Tb = [10 0 -10 10];
physical.qb = [1 -1 -1 1];
%it's heat flux per unit length, we thus need to rescale here
% physical.qb = physical.qb.*[dFine.lElX dFine.lElY dFine.lElX dFine.lElY];
