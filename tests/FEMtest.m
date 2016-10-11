%Test for 2d FEM code
clear all;
addpath('~/matlab/projects/cgrom2d/params')
addpath('~/matlab/projects/cgrom2d/heatFEM')
addpath('~/matlab/projects/cgrom2d/plot')

patchTest = true;
if(patchTest)
    %Test temperature field given by function handle T. FEM solver should lead to the same solution.
    %ONLY MODIFY COEFFICIENTS a, DO NOT MODIFY FUNCTIONAL FORM OF T!!! Otherwise the test will fail
    a = [-1 2 3 -4];
    T = @(x) a(1) + a(2)*x(1) + a(3)*x(2) + a(4)*x(1)*x(2);
    gradT = @(x) [a(2) + a(4)*x(2); a(3) + a(4)*x(1)];
    
    nc = 10;
    domainc = Domain(nc, nc, 1, 1);
    %specify boundary conditions here
    l = 1/nc;
    boundaryCoordinates = [0:l:1, ones(1, nc), (1 - l):(-l):0, zeros(1, nc - 1);...
        zeros(1, nc + 1), l:l:1, ones(1, nc), (1 - l):(-l):l];
    %heat conductivity tensor for each element
    Dc = zeros(2, 2, domainc.nEl);
    for j = 1:domainc.nEl
        %Test is only valid for constant D in the whole domain!
        Dc(:,:,j) = eye(2); %only isotropic material
    end
    for i = 1:4*nc
        physicalc.Tb(i) = T(boundaryCoordinates(:, i));
        qbtemp = - .25*Dc(:, :, 1)*gradT(boundaryCoordinates(:, i));
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
    physicalc.boundaryType = true(1, 4*nc);         %true for essential node, false for natural node
    physicalc.boundaryType(2:nc) = false;           %lower boundary is natural
    physicalc.essentialNodes = domainc.boundaryNodes(physicalc.boundaryType);
    physicalc.naturalNodes = domainc.boundaryNodes(~physicalc.boundaryType);
    domainc = setNodalCoordinates(domainc, physicalc);
    domainc = setBvec(domainc);
    control.plt = false;
    %Assign heat source field
    physicalc.heatSourceField = zeros(domainc.nEl, 1);
    %Force contributions due to heat flux and source
    physicalc.fs = get_heat_source(physicalc.heatSourceField, domainc);
    physicalc.fh = get_flux_force(domainc, physicalc);
    out = heat2d(domainc, physicalc, control, Dc);
    
    for i = 1:domainc.nNodes
        Tcheck(mod(i - 1, nc + 1) + 1, floor((i - 1)/(nc + 1)) + 1) = T(domainc.nodalCoordinates(1:2, i));
    end
    
    testTemperatureField = Tcheck'
    FEMtemperatureField = out.Tff
    figure
    subplot(1, 2, 1)
    pcolor(testTemperatureField);
    title('true temperature field')
    axis square
    subplot(1, 2, 2);
    pcolor(FEMtemperatureField);
    title('FEM temperature field')
    axis square
    if(sqrt(sum(sum((testTemperatureField - FEMtemperatureField).^2)))/numel(testTemperatureField) > 1e-10)
        warning('Patch test for FEM failed')
        difference = sqrt(sum(sum((testTemperatureField - FEMtemperatureField).^2)))/numel(testTemperatureField)
    else
        difference = sqrt(sum(sum((testTemperatureField - FEMtemperatureField).^2)))/numel(testTemperatureField)
        disp('Patch test successful!')
    end
end


convergenceTest = false;
if(convergenceTest)
    % If no parallel pool exists, create one
    N_Threads = 16;
    if isempty(gcp('nocreate'))
        % Create with N_Threads workers
        parpool('local',N_Threads);
    end
    a = [-1 2 3 -4];
    c = 1; %c > 0
    d = 2;
    T = @(x) d*log(norm(x + c)) + a(1) + a(2)*x(1)^2 + a(3)*x(2) + a(4)*x(1)*x(2);
    gradT = @(x) [d*(x(1) + c)/norm(x + c)^2; d*(x(2) + c)/norm(x + c)^2]...
        + [2*a(2)*x(1) + a(4)*x(2); a(3) + a(4)*x(1)];
    
    control.plt = false;
    nSimulations = 30;
    incrementFactor = 3;
    tic;
    for k = 1:nSimulations
        nc = incrementFactor*k;
        domain{k} = Domain(nc, nc, 1, 1);
        %specify boundary conditions here; only essential for this test
        l = 1/nc;
        boundaryCoordinates = [0:l:1, ones(1, nc), (1 - l):(-l):0, zeros(1, nc - 1);...
            zeros(1, nc + 1), l:l:1, ones(1, nc), (1 - l):(-l):l];
        %heat conductivity tensor for each element
        Dc = zeros(2, 2, domain{k}.nEl);
        for j = 1:domain{k}.nEl
            %Test is only valid for constant D in the whole domain!
            Dc(:,:,j) = eye(2); %only isotropic material
        end
        for i = 1:4*nc
            physical{k}.Tb(i) = T(boundaryCoordinates(:, i));
            qbtemp = - .25*Dc(:, :, 1)*gradT(boundaryCoordinates(:, i));
            %projection along normal vectors of domain boundaries
            if i <= nc
                %bottom
                physical{k}.qb(i) = qbtemp(2);
            elseif(mod(i, nc + 1) == 0 && i < (nc + 1)^2)
                %right
                physical{k}.qb(i) = -qbtemp(1);
            elseif(i > nc*(nc + 1))
                %top
                physical{k}.qb(i) = -qbtemp(2);
            elseif(mod(i, nc + 1) == 1 && i > 1)
                %left
                physical{k}.qb(i) = qbtemp(1);
            end
        end
        physical{k}.boundaryType = true(1, 4*nc);         %true for essential node, false for natural node
        %     physical{k}.boundaryType(2:nc) = false;           %lower boundary is natural
        physical{k}.essentialNodes = domain{k}.boundaryNodes(physical{k}.boundaryType);
        physical{k}.naturalNodes = domain{k}.boundaryNodes(~physical{k}.boundaryType);
        domain{k} = setNodalCoordinates(domain{k}, physical{k});
        domain{k} = setBvec(domain{k});
        %Assign heat source field
        physical{k}.heatSourceField = zeros(domain{k}.nEl, 1);
        %Force contributions due to heat flux and source
        physical{k}.fs = get_heat_source(physical{k}.heatSourceField, domain{k});
        physical{k}.fh = get_flux_force(domain{k}, physical{k});
        for i = 1:domain{k}.nNodes
            Tcheck{k}(mod(i - 1, nc + 1) + 1, floor((i - 1)/(nc + 1)) + 1) = T(domain{k}.nodalCoordinates(1:2, i));
        end
        testTemperatureField{k} = Tcheck{k}';
    end
    t1 = toc;
    parfor k = 1:nSimulations
    out = heat2d(domain{k}, physical{k}, control, Dc);
        FEMtemperatureField{k} = out.Tff;
        difference(k) = sqrt(sum(sum((testTemperatureField{k} - FEMtemperatureField{k}).^2)))/numel(testTemperatureField{k});
        nElementsX(k) = domain{k}.nElX;
    end
    figure
    loglog(nElementsX, difference)
    xlabel('Number of elements')
    ylabel('Root mean square difference')
    title('Convergence to true solution')
    
    
end