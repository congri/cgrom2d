%performance tests

a = rand(size(theta_cf.W, 1), 1);
nRuns = 10000;
tic
for i = 1:nRuns
    b = theta_cf.W'*a;
end
t1 = toc

tic
WT = theta_cf.W';
for i = 1:nRuns
    b = WT*a;
end
t2 = toc