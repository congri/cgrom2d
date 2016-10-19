function [Out] = heat2d(domain, D)
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
Out.globalForce = get_glob_force(domain, Out.localStiffness);

%Finally solving the equation system
Out.naturalTemperatures = Out.globalStiffness\Out.globalForce;

%Temperature field
Tf = zeros(domain.nNodes, 1);
Tf(domain.id) = Out.naturalTemperatures;
Tff = zeros(domain.nElX + 1, domain.nElY + 1);

for i = 1:domain.nNodes
    Tff(i) = Tf(i);
    if(any(i == domain.essentialNodes))
        %node i is essential
        Tff(i) = domain.essentialTemperatures(i);
    end
end
Tff = Tff';
Out.Tff = Tff;


end