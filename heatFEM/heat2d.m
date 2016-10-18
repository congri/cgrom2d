function [Out] = heat2d(domain, physical, control, D)
%2D heat conduction main function
%Gives back temperature on point x

%get_loc_stiff as nested function for performance 
function [k] = get_loc_stiff2(Bvec, D)
    %Gives the local stiffness matrix

    Z = zeros(2, 2);
    Dmat = [D Z Z Z; Z D Z Z; Z Z D Z; Z Z Z D];

    k = Bvec'*Dmat*Bvec;
end


%Compute local stiffness matrices, once and for all
Out.localStiffness = zeros(4, 4, domain.nEl);
for e = 1:domain.nEl
    Out.localStiffness(:, :, e) = get_loc_stiff2(domain.Bvec(:, :, e), D(:, :, e));
end

%Global stiffness matrix
Out.globalStiffness = get_glob_stiff2(domain, Out.localStiffness);
%Global force vector
Out.globalForce = get_glob_force(domain, physical, Out.localStiffness);

%Finally solving the equation system
Out.naturalTemperatures = Out.globalStiffness\Out.globalForce;

%Temperature field
Tf = zeros(domain.totalNodeNumber,1);

Tf(domain.id) = Out.naturalTemperatures;
% if(strcmp(domain.boundaries(1), 'essential'))
%     Tf(1:(domain.nElX + 1)) = physical.Tb(1);
% end
% if(strcmp(domain.boundaries(2),'essential'))
%    Tf((domain.nElX + 1):(domain.nElX + 1):(domain.nElX + 1)*(domain.nElY + 1)) = physical.Tb(2);
% end
% if(strcmp(domain.boundaries(3),'essential'))
%     Tf((domain.nElY*(domain.nElX + 1) + 1):(domain.nElY + 1)*(domain.nElX + 1)) = physical.Tb(3);
% end
% if(strcmp(domain.boundaries(4),'essential'))
%    Tf(1:(domain.nElX + 1):(domain.nElX + 1)*(domain.nElY + 1)) = physical.Tb(4); 
% end

Tff = zeros(domain.nElX + 1, domain.nElY + 1);

for i = 1:domain.totalNodeNumber
   Tff(i) = Tf(i);
   if(~isnan(domain.nodalCoordinates(4, i)))
       Tff(i) = domain.nodalCoordinates(4, i);
   end
end
Tff = Tff';
Out.Tff = Tff;

%Find temperature Tx on input point x
[X, Y] = meshgrid(linspace(0, domain.lx, domain.nElX + 1), linspace(0, domain.ly, domain.nElY + 1));

%Plot?
if(control.plt)
    plotHeatMap(X, Y, Tff, domain, physical);
    %[Qx, Qy] = plotflux(X, Y, Tff, domain, D);
end
end