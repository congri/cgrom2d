function [s] = specificSurface(lambda, phase, fineData, nElc, nElf)
%The specific surface is given by the derivative of the 2-point correlation function at 0, see
%Torquato 2.2.5, in particular equations 2.32 and 2.36.

S0 = twoPointCorrelation(lambda, 0, 'x', phase, fineData, nElc, nElf);
S1 = .5*twoPointCorrelation(lambda, 1, 'x', phase, fineData, nElc, nElf) +...
        .5*twoPointCorrelation(lambda, 1, 'y', phase, fineData, nElc, nElf);

%pixel edge length
%Wrong if domain length ~=1!!!
l_pixel_x = 1/nElf(1);
l_pixel_y = 1/nElf(2);

dx = (S1 - S0)/l_pixel_x;
dy = (S1 - S0)/l_pixel_y;

%is this correct?
dr = .5*dx + .5*dy;
s = -pi*dr;     %Torquato 2.36 for 2D

end

