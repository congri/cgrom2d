function [out] = linPathParams(cond, maxPathLength, fineData, domainc, domainf, param)
%Gives back the parameters a and b from the theoretical lineal path model
%L(z) = a*exp(b*z)

L = zeros(maxPathLength + 1, 1);
nElc = [domainc.nElX domainc.nElY];
nElf = [domainf.nElX domainf.nElY];
for i = 0:maxPathLength
    L(i + 1) = .5*linealPath(cond, i, 'x', 1, fineData, nElc, nElf)...
        + .5*linealPath(cond, i, 'y', 1, fineData, nElc, nElf);
end

d = (0:maxPathLength)';
f = fit(d, L, 'exp1');
if strcmp(param, 'a')
    out = f.a;
elseif strcmp(param, 'b')
    out = f.b;
else
    error('Which Lineal path parameter?')
end

