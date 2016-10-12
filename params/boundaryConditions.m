%generate boundary conditions here

%% Temperature field and gradient generating the boundary conditions
a = [-4 5 7 -4];
c = .3; %c > 0
d = 4;
T = @(x) d*log(norm(x + c)) + a(1) + a(2)*x(1)^2 + a(3)*x(2) + a(4)*x(1)*x(2);
gradT = @(x) [d*(x(1) + c)/norm(x + c)^2; d*(x(2) + c)/norm(x + c)^2]...
    + [2*a(2)*x(1) + a(4)*x(2); a(3) + a(4)*x(1)];


%% coarse boundary values
lc = 1/nc;
boundaryCoordinatesCoarse = [0:lc:1, ones(1, nc), (1 - lc):(-lc):0, zeros(1, nc - 1);...
    zeros(1, nc + 1), lc:lc:1, ones(1, nc), (1 - lc):(-lc):lc];
for i = 1:4*nc
    physicalc.Tb(i) = T(boundaryCoordinatesCoarse(:, i));
    qbtemp = - .25*eye(2)*gradT(boundaryCoordinatesCoarse(:, i)); %the factor .25 is not explicable,
    %it must be missing somwhere in the FEM code
    
    %projection along normal vectors of domain boundaries
    if i <= nc
        %bottom
        physicalc.qb(i) = qbtemp(2);
    elseif(mod(i, nc + 1) == 0 && i < (nc + 1)^2)
        %right
        physicalc.qb(i) = -qbtemp(1);
    elseif(i > nc*(nc + 1))
        %top
        physicalc.qb(i) = -qbtemp(2);
    elseif(mod(i, nc + 1) == 1 && i > 1)
        %left
        physicalc.qb(i) = qbtemp(1);
    end
end
%ATTENTION: BOUNDARY TYPES OF COARSE AND FINE MODELS HAVE TO BE THE SAME ON CORRESPONDING BOUNDARIES
physicalc.boundaryType = true(1, 4*nc);         %true for essential node, false for natural node
physicalc.boundaryType(2:nc) = false;           %lower boundary is natural
physicalc.boundaryType((nc + 2):(2*nc)) = false;%right boundary is natural
physicalc.essentialNodes = domainc.boundaryNodes(physicalc.boundaryType);
physicalc.naturalNodes = domainc.boundaryNodes(~physicalc.boundaryType);
%Assign heat source field
physicalc.heatSourceField = zeros(domainc.nEl, 1);
%nodal coordinates must be set here
domainc = setNodalCoordinates(domainc, physicalc);
domainc = setBvec(domainc);
%Force contributions due to heat flux and source
physicalc.fs = get_heat_source(physicalc.heatSourceField, domainc);
physicalc.fh = get_flux_force(domainc, physicalc);
physicalc.floc = physicalc.fs + physicalc.fh;


%% file boundary values
lf = 1/nf;
boundaryCoordinatesFine = [0:lf:1, ones(1, nf), (1 - lf):(-lf):0, zeros(1, nf - 1);...
    zeros(1, nf + 1), lf:lf:1, ones(1, nf), (1 - lf):(-lf):lf];
for i = 1:4*nf
    physicalf.Tb(i) = T(boundaryCoordinatesFine(:, i));
    qbtemp = - .25*eye(2)*gradT(boundaryCoordinatesFine(:, i)); %the factor .25 is not explicable,
    %it must be missing somwhere in the FEM code
    
    %projection along normal vectors of domain boundaries
    if i <= nf
        %bottom
        physicalf.qb(i) = qbtemp(2);
    elseif(mod(i, nf + 1) == 0 && i < (nf + 1)^2)
        %right
        physicalf.qb(i) = -qbtemp(1);
    elseif(i > nf*(nf + 1))
        %top
        physicalf.qb(i) = -qbtemp(2);
    elseif(mod(i, nf + 1) == 1 && i > 1)
        %left
        physicalf.qb(i) = qbtemp(1);
    end
end
%ATTENTION: BOUNDARY TYPES OF COARSE AND FINE MODELS HAVE TO BE THE SAME ON CORRESPONDING BOUNDARIES
physicalf.boundaryType = true(1, 4*nf);         %true for essential node, false for natural node
physicalf.boundaryType(2:nf) = false;           %lower boundary is natural
physicalf.boundaryType((nf + 2):(2*nf)) = false;%right boundary is natural
physicalf.essentialNodes = domainf.boundaryNodes(physicalf.boundaryType);
physicalf.naturalNodes = domainf.boundaryNodes(~physicalf.boundaryType);
%Assign heat source field
physicalf.heatSourceField = zeros(domainf.nEl, 1);
%nodal coordinates must be set here
domainf = setNodalCoordinates(domainf, physicalf);
domainf = setBvec(domainf);
%Force contributions due to heat flux and source
physicalf.fs = get_heat_source(physicalf.heatSourceField, domainf);
physicalf.fh = get_flux_force(domainf, physicalf);
physicalf.floc = physicalf.fs + physicalf.fh;