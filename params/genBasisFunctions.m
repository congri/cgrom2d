%% Some predefined basis functions for linear model p_c

% phi = linPath;
% phi{1} = lpa;
% phi{end + 1} = lpb;
% phi{end + 1} = mps;
% phi{end + 1} = vps;
% phi{end + 1} = phi_1;
% phi{end + 1} = phi_2;
% phi{end + 1} = phi_3;

%harmonic, arhithmetic and geometric exponentiated
% minexp = 0;
% maxexp = 5;
% expincr = .5;
% exponents = minexp:expincr:maxexp;
% nExponents = length(exponents);
% i = 1;
% for ex = exponents
%     phi{i} = @(lambda) phi_2(lambda)^ex;
%     phi{i + nExponents} = @(lambda) phi_1(lambda)^ex;
%     phi{i + 2*nExponents} = @(lambda) phi_3(lambda)^ex;
%     i = i + 1;
% end

%generalized means
min_z = -1;
max_z = 1;
z_incr = .1;
z = min_z:z_incr:max_z;
i = 1;
for zz = z
    phi{i} = @(lambda) generalizedMean(lambda, zz);
    i = i + 1;
end
%log of generalized means
for zz = z
    phi{i} = @(lambda) log(generalizedMean(lambda, zz));
    i = i + 1;
end

dLinPathMax = 80;
dLinPathMin = 0;
dLinPathIncr = 4;
nElc = [domainc.nElX domainc.nElY];
nElf = [domainf.nElX domainf.nElY];
% for d = dLinPathMin:dLinPathIncr:dLinPathMax
%     phi{end + 1} = @(lambda) .5*linealPath(lambda, d, 'x', 2, fineData, nElc, nElf) +...
%     .5*linealPath(lambda, d, 'y', 2, fineData, nElc, nElf);
% end
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    phi{end + 1} = @(lambda) linealPath(lambda, d, 'x', 2, fineData, nElc, nElf);
end
for d = dLinPathMin:dLinPathIncr:dLinPathMax
    phi{end + 1} = @(lambda) linealPath(lambda, d, 'y', 2, fineData, nElc, nElf);
end

d2pointCorrMax = 100;
d2pointCorrMin = 2;
d2pointCorrIncr = 1;
i = 1;
% for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
%     phi{end + 1} = @(lambda) .5*twoPointCorrelation(lambda, d, 'x', 2, fineData, nElc, nElf) +...
%         .5*twoPointCorrelation(lambda, d, 'y', 2, fineData, nElc, nElf);
%     i = i + 1;
% end
for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
    %distinct in x and y direction
    phi{end + 1} = @(lambda) twoPointCorrelation(lambda, d, 'x', 2, fineData, nElc, nElf);
end
for d = d2pointCorrMin:d2pointCorrIncr:d2pointCorrMax
    %distinct in x and y direction
    phi{end + 1} = @(lambda) twoPointCorrelation(lambda, d, 'y', 2, fineData, nElc, nElf);
end

pathLengths = (0:4:60)';
% lpa = @(lambda) linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'a');
lpb = @(lambda) abs(linPathParams(lambda, pathLengths, fineData, domainc, domainf, 'b'));

mps = @(lambda) meanPoreSize(lambda, 2, fineData, nElc, nElf, 'mean')/(nElf(1)/2);
vps = @(lambda) sqrt(meanPoreSize(lambda, 2, fineData, nElc, nElf, 'var'))/(nElf(1)/2);

phi{end + 1} = mps;
phi{end + 1} = vps;
phi{end + 1} = @(lambda) specificSurface(lambda, 2, fineData, nElc, nElf);
phi{end + 1} = lpb;
%This gives mean euclidean distance to nearest high conducting pixel
phi{end + 1} = @(lambda) mean(mean(bwdist(reshape(lambda > fineData.lo, sqrt(length(lambda)), sqrt(length(lambda))))));
%This gives mean euclidean distance to nearest low conducting pixel
phi{end + 1} = @(lambda) mean(mean(bwdist(reshape(lambda < fineData.up, sqrt(length(lambda)), sqrt(length(lambda))))));
%Counts the number of separate high conducting blobs
phi{end + 1} = @(lambda) numberOfObjects(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))), fineData, 'hi');
%Counts the number of separate low conducting blobs
phi{end + 1} = @(lambda) numberOfObjects(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))), fineData, 'lo');

%% High conducting phase
%means
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Area', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'ConvexArea', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Eccentricity', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Extent', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MajorAxisLength', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MinorAxisLength', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Orientation', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Perimeter', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Solidity', 'mean');
%Variances
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Area', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'ConvexArea', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Eccentricity', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Extent', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MajorAxisLength', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MinorAxisLength', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Orientation', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Perimeter', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Solidity', 'var');
%Maxima
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Area', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'ConvexArea', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Eccentricity', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Extent', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MajorAxisLength', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MinorAxisLength', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Orientation', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Perimeter', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Solidity', 'max');
%Minima
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Area', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'ConvexArea', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Eccentricity', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Extent', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MajorAxisLength', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MinorAxisLength', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Orientation', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Perimeter', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Solidity', 'min');
%standard deviations
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Area', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'ConvexArea', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Eccentricity', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Extent', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MajorAxisLength', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'MinorAxisLength', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Orientation', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Perimeter', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'Solidity', 'var').^.5;
%maximum bubble extents
phi{end + 1} = @(lambda) maxExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'x');
phi{end + 1} = @(lambda) maxExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'y');
%mean bubble extents
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'x', 'mean');
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'y', 'mean');
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'x', 'var');
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'hi', 'y', 'var');

%% Low conducting phase
%means
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Area', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'ConvexArea', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Eccentricity', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Extent', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MajorAxisLength', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MinorAxisLength', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Orientation', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Perimeter', 'mean');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Solidity', 'mean');
%Variances
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Area', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'ConvexArea', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Eccentricity', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Extent', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MajorAxisLength', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MinorAxisLength', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Orientation', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Perimeter', 'var');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Solidity', 'var');
%Maxima
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Area', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'ConvexArea', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Eccentricity', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Extent', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MajorAxisLength', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MinorAxisLength', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Orientation', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Perimeter', 'max');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Solidity', 'max');
%Minima
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Area', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'ConvexArea', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Eccentricity', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Extent', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MajorAxisLength', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MinorAxisLength', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Orientation', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Perimeter', 'min');
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Solidity', 'min');
%standard deviations
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Area', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'ConvexArea', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Eccentricity', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Extent', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MajorAxisLength', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'MinorAxisLength', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Orientation', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Perimeter', 'var').^.5;
phi{end + 1} = @(lambda) meanImageProps(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'Solidity', 'var').^.5;
%maximum bubble extents
phi{end + 1} = @(lambda) maxExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'x');
phi{end + 1} = @(lambda) maxExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'y');
%mean bubble extents
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'x', 'mean');
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'y', 'mean');
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'x', 'var');
phi{end + 1} = @(lambda) meanExtent(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))),...
    fineData, 'lo', 'y', 'var');



%Energy of conductivity field if it were ising model
phi{end + 1} = @(lambda) isingEnergy(reshape(lambda, sqrt(length(lambda)), sqrt(length(lambda))));
%constant
phi{end + 1} = @(lambda) 1;

%max/min
phi{end + 1} = @(lambda) max(lambda);
phi{end + 1} = @(lambda) min(lambda);



%% Sum of random pixels
%Fine elements per coarse element in x and y directions
% nRandPixelFeatures = 100;
% xc = nElf(1)/nElc(1);
% yc = nElf(2)/nElc(2);
% finePerCoarse = xc*yc;
% randomPixels = randi([0 1], nRandPixelFeatures, finePerCoarse);
% for i = 1:nRandPixelFeatures
%     phi{end + 1} = @(lambda) randomPixels(i, :)*lambda;
% end




%Pixel grid
%Fine elements per coarse element in x and y directions
% xc = nElf(1)/nElc(1);
% yc = nElf(2)/nElc(2);
% firstRow = 33;
% lastRow = yc - 32;
% firstCol = 33;
% lastCol = xc - 32;
% pixelIncr = 64;
% for row = firstRow:pixelIncr:lastRow
%     for col = firstCol:pixelIncr:lastCol
%         phi{end + 1} = @(lambda) lambda(col + xc*(row - 1));
%     end
% end
% 
nBasis = numel(phi);


%print feature function strings to file
delete('./data/basisFunctions.txt');
fid = fopen('./data/basisFunctions.txt','wt');
for i = 1:length(phi)
    phi_str = func2str(phi{i});
    phi_str = strcat(phi_str, '\n');
    fprintf(fid, phi_str);
end
fclose(fid);



