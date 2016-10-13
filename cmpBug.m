%Colormap bug?
%only bugged on cluster

nData = 10;
A = rand(nData);
cmp = [0 0 .5; 0 0 .7; 0 0 1];
[X, Y] = meshgrid(linspace(1, nData, nData));
figure
[~, cf] = contourf(X, Y, A)
cf.LineWidth = 1;
colormap(cmp)
colorbar
