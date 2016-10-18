function [nc] = get_coord(domain, physical)
%Gives nodal coordinates in the first two rows and equation number from
%global node number in the third row. Temperature of essential boundaries
%is given in the fourth row, heat flux on natural boundaries in the fifth
%Assign global node coordinates and equation numbers
%In clockwise direction, the first node of every side is considered to belong to the boundary. The
%last node is considered to belong to the next boundary. E.g. on a grid 5x5 nodes, nodes 1 - 4
%belong to the lower boundary, nodes 5, 10, 15, 20 to the right, nodes 25, 24, 23, 22 to the upper
%and 21, 16, 11, 6 to the left boundary.

x = 0;
y = 0;
j = 1;  %equation number index
k = 1;  %boundary node index
nc = NaN*zeros(4, domain.totalNodeNumber);
for i = 1:domain.totalNodeNumber
    nc(1, i) = x;
    nc(2, i) = y;

    if(any(physical.essentialNodes == i))
        %essential node, no equation number assigned
        nc(3, i) = 0;
    else
        %Assign equation number j
        nc(3, i) = j;
        j = j + 1;
    end
    
    x = x + domain.lElX;
    %reset x to 0 on last node of each row; increase y
    if mod(i, domain.nElX + 1) == 0
        x = 0;
        y = y + domain.lElY;
    end
end
nc(4, domain.essentialNodes) = physical.Tb(domain.boundaryType);

end

