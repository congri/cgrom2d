function [Tf_mean_mat, Tf_std_mat] = plotResult(theta_c, theta_cf, domainc, domainf, phi)
%We sample the resulting distribution p(y|x, theta_c, theta_cf) here and compare to the true
%solution

tic;
%load finescale data from last optimization
load('./data/fineData/testData')


%design Matrix for p_c
Phi = designMatrix(phi, condTest, domainf, domainc);

%% Sample from p_c
nSamples_p_c = 100;
Xsamples = mvnrnd(Phi*theta_c.theta, (theta_c.sigma^2)*eye(domainc.nEl), nSamples_p_c)';
LambdaSamples = exp(Xsamples);

%% Run coarse model and sample from p_cf
Tc = zeros(domainc.nNodes, nSamples_p_c);
Tfinterp = zeros(domainf.nNodes, nSamples_p_c);
for i = 1:nSamples_p_c
    D = zeros(2, 2, domainc.nEl);
    for e = 1:domainc.nEl
        D(:, :, e) = LambdaSamples(e, i)*eye(2);
    end
    FEMout = heat2d(domainc, D);
    Tctemp = FEMout.Tff';
    Tc(:, i) = Tctemp(:);
    
    %sample from p_cf
    mu_cf = theta_cf.mu + theta_cf.W*Tc(:, i);
    Tfinterp(:, i) = mvnrnd(mu_cf', theta_cf.S);
end

%% Plot
Tf_mean = mean(Tfinterp, 2);
Tf_mean_mat = reshape(Tf_mean, domainf.nElX + 1, domainf.nElY + 1);
Tf_mean_mat = Tf_mean_mat';
Tf_std = std(Tfinterp');
Tf_std_mat = reshape(Tf_std, domainf.nElX + 1, domainf.nElY + 1);
Tf_std_mat = Tf_std_mat';
[Xcoord, Ycoord] = meshgrid(linspace(0, 1, domainf.nElX + 1), linspace(0, 1, domainf.nElY + 1));

Tf_true_mat = reshape(TfTest, domainf.nElX + 1, domainf.nElY + 1);
Tf_true_mat = Tf_true_mat';

LambdacMean = mean(LambdaSamples, 2)
LambdacMean = reshape(LambdacMean, domainc.nElX, domainc.nElY)';
LambdacMeanPlot = zeros(size(LambdacMean) + 1);
LambdacMeanPlot(1:(end - 1), 1:(end - 1)) = LambdacMean;
[LambdacX, LambdacY] = meshgrid(linspace(0, 1, domainc.nElX + 1), linspace(0, 1, domainc.nElY + 1));

Lambdaf = reshape(condTest, domainf.nElX, domainf.nElY)';
[LambdafX, LambdafY] = meshgrid(linspace(0, 1, domainf.nElX + 1), linspace(0, 1, domainf.nElY + 1));
LambdafPlot = zeros(size(Lambdaf) + 1);
LambdafPlot(1:(end - 1), 1:(end - 1)) = Lambdaf;

%Error measure
err = sqrt((1./(2*Tf_std'.^2)).*(TfTest - Tf_mean).^2);
err_mat = reshape(err, domainf.nElX + 1, domainf.nElY + 1);
err_mat = err_mat';


cmin = min([min(min(Tf_mean_mat)), min(min(Tf_true_mat))]);
cmax = max([max(max(Tf_mean_mat)), max(max(Tf_true_mat))]);

cmp = inferno();
figure
subplot(3,2,2)
contourf(Xcoord, Ycoord, Tf_mean_mat, 256, 'linestyle', 'none')
xlabel('x')
ylabel('y')
title('Pred. mean')
axis square
colormap(cmp)
colorbar
hold
plot([0 1], [.5 .5], 'w', 'linewidth', 3)
% caxis([cmin, cmax])

subplot(3,2,1)
contourf(Xcoord, Ycoord, Tf_true_mat, 256, 'linestyle', 'none')
xlabel('x')
ylabel('y')
title('True finescale output')
axis square
colormap(cmp)
colorbar
hold
plot([0 1], [.5 .5], 'w', 'linewidth', 3)
% caxis([cmin, cmax])

subplot(3,2,3)
pcolor(LambdacX, LambdacY, LambdacMeanPlot)
title('Mean eff. conductivity')
xlabel('x')
ylabel('y')
grid off
colormap(cmp)
colorbar
caxis([min(min(LambdacMeanPlot(1:(end - 1), 1:(end - 1)))) max(max(LambdacMeanPlot(1:(end - 1), 1:(end - 1))))])
axis square

subplot(3,2,4)
contourf(Xcoord, Ycoord, Tf_std_mat, 256, 'linestyle', 'none')
title('Pred. standard deviation')
xlabel('x')
ylabel('y')
axis square
colormap(cmp)
colorbar

subplot(3,2,5)
pcolor(LambdafX, LambdafY, LambdafPlot)
title('True conductivity')
xlabel('x')
ylabel('y')
grid off
colormap(cmp)
colorbar
caxis([min(min(LambdafPlot(1:(end - 1), 1:(end - 1)))), max(max(LambdafPlot(1:(end - 1), 1:(end - 1))))])
axis square

% subplot(3,2,6)
% surf(Xcoord, Ycoord, Tf_mean_mat)
% hold on
% surf(Xcoord, Ycoord, Tf_true_mat)
% axis square

subplot(3,2,6)
contourf(Xcoord, Ycoord, err_mat, 256, 'linestyle', 'none')
title('Error measure')
xlabel('x')
ylabel('y')
colorbar
axis square

%line cut through middle
s = size(Tf_true_mat);
T_true_cut = Tf_true_mat(floor(s(1)/2), :);
T_pred_cut = Tf_mean_mat(floor(s(1)/2), :);
pred_std = Tf_std_mat(floor(s(1)/2), :);
figure
x = linspace(0, 1, s(2));
shadedErrorBar(x, T_pred_cut, pred_std);
hold
plot(x, T_true_cut, 'linewidth', 3)
axis square
xlabel('x')
ylabel('T')

runtime = toc
end






