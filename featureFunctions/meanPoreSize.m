function [m] = meanPoreSize(lambda, phase, fineData, nElc, nElf)
%Computes mean pore size density
[p_delta_x, p_delta] = poreSizeDensity(lambda, phase, fineData, nElc, nElf);

m = p_delta_x'*p_delta;

end

