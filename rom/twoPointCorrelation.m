function [p] = twoPointCorrelation(lambdak, distance, dir, phase, fineData, domainc, domainf)
%two-point correlation function of materials
%   lambdak:        fine conductivities in coarse element k
%   distance:     length of distance in units of fine elements (pixels). Must be smaller than sqrt(nFine/nCoarse)
%   dir:            x or y direction
%   phase:          Material phase, 1 or 2 for binary
%   fineData:       distribution type of fine conductivities. Must be binary!
%
%   Output:
%       p:          2-point correlation function to given parameters, see Torquato p. 25

%Check
% assert(strcmp(fineData.dist, 'binary'), 'Error: 2-point correlation is only a possible basis function if conductivities are binary')

phaseConductivity = [fineData.lo; fineData.up];    %conductivities of the two phases

%Fine elements per coarse element in x and y directions
xc = domainf.nElX/domainc.nElX;
yc = domainf.nElY/domainc.nElY;

if dir == 'x'
    %Maximal number of line segments of length distance that can be dropped in a single row/column in
    %x/y direction
    maxX = xc - distance;
    assert(maxX >= 1, 'Error: maximal number of possible line segments per row must be at least 1')
    samePhase = 0;
    i = 1;
    for row = 1:yc
        for col = 1:maxX
            if(lambdak(i) == phaseConductivity(phase) && lambdak(i + distance) == phaseConductivity(phase))
                samePhase = samePhase + 1;
            end
            i = i + 1;
        end
        i = i + distance;   %jump over to first column in next row
    end
    p = samePhase/(maxX*yc);
    
elseif dir == 'y'
    %Maximal number of line segments of length pathLength that can be dropped in a single row/column in
    %x/y direction
    maxY = yc - distance;
    assert(maxY >= 1, 'Error: maximal number of possible line segments per column must be at least 1')
    samePhase = 0;
    i = 1;
    for row = 1:maxY
        for col = 1:xc
            if(lambdak(i) == phaseConductivity(phase) && lambdak(i + distance*xc) == phaseConductivity(phase))
               samePhase = samePhase + 1; 
            end
            i = i + 1;
        end
        %no jump needed here?
    end
    p = samePhase/(maxY*xc);
    
else
    error('Unknown direction for twoPointCorrelation function')
end



end

