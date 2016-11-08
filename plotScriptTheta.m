%Plotscript to plot theta params
dLinPathMax = 160;
dLinPathMin = 0;
dLinPathIncr = 8;
d2pointCorrMax = 160;
d2pointCorrMin = 8;
d2pointCorrIncr = 8;
subplot(3, 2, 1)
pl = plot(dLinPathMin:dLinPathIncr:dLinPathMax,...
    theta_c.theta(1:(((dLinPathMax - dLinPathMin)/dLinPathIncr) + 1)), 'linewidth', 2);
hold on
p2p = plot(d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax,...
    theta_c.theta((((dLinPathMax - dLinPathMin)/dLinPathIncr) + 2):end), 'linewidth', 2);
xlim([min([dLinPathMin, d2pointCorrMin]) max([dLinPathMax d2pointCorrMax])])
legend([pl, p2p], 'lineal path', '2-point')
xlabel('pixels distance')
ylabel('\theta_i')