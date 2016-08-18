%main script for 2d coarse-graining

addpath('./params')
addpath('./heatFEM')
addpath('./plot')
addpath('./rom')

clear all;
%load params
params;

%genereate full order model data
[x, PhiArray, Tf] = genFineData(dFine, dCoarse, fineCond, phi, physF);


% control.plt = 0;
% lambdaC = mean(mean(lambdaUnit))*ones(dCoarse.nEl, 1);
% Dcoarse = zeros(2,2,dCoarse.nEl);
% %coarse model
% for i = 1:dCoarse.nEl
%     Dcoarse(:,:,i) = lambdaC(i)*eye(2); %only isotropic material
% end
% [coarseOut] = heat2d(dCoarse, physC, control, Dcoarse)
% 
% %fine model
% lambdaF = repmat(lambdaUnit, nF/2, nF/2);
% lambdaF = lambdaF(:);
% 
% Dfine = zeros(2,2,dFine.nEl);
% for i = 1:dFine.nEl
%     Dfine(:,:,i) = lambdaF(i)*eye(2); %only isotropic material
% end
% [fineOut] = heat2d(dFine, physF, control, Dfine)

%store handle to every q_i in a cell array lq
lq = cell(fineCond.nSamples, 1);











f = figure('name','2D heat map');
set(f, 'Position', [10, 350, 720, 540]);
subplot(2,2,1)
[Xf, Yf] = meshgrid(linspace(0, dFine.lx, dFine.nElX + 1), linspace(0, dFine.ly, dFine.nElY + 1));
contourf(Xf, Yf, Tf(:,:,1), 256,'LineColor','none');
axis square
title('2D heat map');
xlabel('x');
ylabel('y');
cb = colorbar;
ylabel(cb,'T');
set(gca,'FontSize',14) 
colormap('jet')
subplot(2,2,2)
[Xf, Yf] = meshgrid(linspace(0, dFine.lx, dFine.nElX + 1), linspace(0, dFine.ly, dFine.nElY + 1));

contourf(Xf, Yf, Tf(:,:,2), 256,'LineColor','none');
axis square
title('2D heat map');
xlabel('x');
ylabel('y');
cb = colorbar;
ylabel(cb,'T');
set(gca,'FontSize',14) 
colormap('jet')
