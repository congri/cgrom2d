%main script for 2d coarse-graining

addpath('./params')
addpath('./heatFEM')
addpath('./plot')

clear all;
%load params
params;

control.plt = 1;

lambda = ones(dFine.nEl, 1);
D = zeros(2,2,dFine.nEl);
for i = 1:dFine.nEl
    D(:,:,i) = lambda(i)*eye(2); %only isotropic material
end
[Out] = heat2d(dFine, physical, control, D)