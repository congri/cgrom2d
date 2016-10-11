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
nc = NaN*zeros(5,domain.totalNodeNumber);
for i = 1:domain.totalNodeNumber
    nc(1, i) = x;
    nc(2, i) = y;
    
%     if(i == 1) %lower left corner
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(1),'essential') || strcmp(domain.boundaries(4),'essential'))
%             %if any of the two boundaries is essential, the node is essential
%             nc(3,i) = 0;
%             if(strcmp(domain.boundaries(1),'essential') && strcmp(domain.boundaries(4),'essential'))
%                 nc(4,i) = 0.5*(physical.Tb(1) + physical.Tb(4)); %We take the mean value of boundaries
%             elseif(strcmp(domain.boundaries(1),'essential'))
%                 nc(4,i) = physical.Tb(1);
%             elseif(strcmp(domain.boundaries(4),'essential'))
%                 nc(4,i) = physical.Tb(4);
%             else
%                 error('No essential boundary?');
%             end
%         elseif(strcmp(domain.boundaries(1),'natural') && strcmp(domain.boundaries(4), 'natural'))
%             %if both boundaries are natural, the node is natural
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = 0.5*(physical.qb(1) + physical.qb(4)); %we take the mean value
%         elseif(isnan(nc(3,i)))
%             error('Neither essential nor natural boundary');
%         end
%     end
%     
%     if(i == domain.nElX + 1) %lower right corner
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(1),'essential') || strcmp(domain.boundaries(2),'essential'))
%             %if any of the two boundaries is essential, the node is essential
%             nc(3,i) = 0;
%             if(strcmp(domain.boundaries(1),'essential') && strcmp(domain.boundaries(2),'essential'))
%                 nc(4,i) = 0.5*(physical.Tb(1) + physical.Tb(2)); %We take the mean value of boundaries
%             elseif(strcmp(domain.boundaries(1),'essential'))
%                 nc(4,i) = physical.Tb(1);
%             elseif(strcmp(domain.boundaries(2),'essential'))
%                 nc(4,i) = physical.Tb(2);
%             else
%                 error('No essential boundary?');
%             end
%         elseif(strcmp(domain.boundaries(1),'natural') && strcmp(domain.boundaries(2), 'natural'))
%             %if both boundaries are natural, the node is natural
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = 0.5*(physical.qb(1) + physical.qb(2)); %we take the mean value
%         elseif(isnan(nc(3,i)))
%             error('Neither essential nor natural boundary');
%         end
%     end
%     
%     if(i == (domain.nElX + 1)*(domain.nElY + 1)) %upper right corner
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(3),'essential') || strcmp(domain.boundaries(2),'essential'))
%             %if any of the two boundaries is essential, the node is essential
%             nc(3,i) = 0;
%             if(strcmp(domain.boundaries(3),'essential') && strcmp(domain.boundaries(2),'essential'))
%                 nc(4,i) = 0.5*(physical.Tb(3) + physical.Tb(2)); %We take the mean value of boundaries
%             elseif(strcmp(domain.boundaries(3),'essential'))
%                 nc(4,i) = physical.Tb(3);
%             elseif(strcmp(domain.boundaries(2),'essential'))
%                 nc(4,i) = physical.Tb(2);
%             else
%                 error('No essential boundary?');
%             end
%         elseif(strcmp(domain.boundaries(3),'natural') && strcmp(domain.boundaries(2), 'natural'))
%             %if both boundaries are natural, the node is natural
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = 0.5*(physical.qb(3) + physical.qb(2)); %we take the mean value
%         elseif(isnan(nc(3,i)))
%             error('Neither essential nor natural boundary');
%         end
%     end
%     
%     if(i == (domain.nElX + 1)*domain.nElY + 1) %upper left corner
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(3),'essential') || strcmp(domain.boundaries(4),'essential'))
%             %if any of the two boundaries is essential, the node is essential
%             nc(3,i) = 0;
%             if(strcmp(domain.boundaries(3),'essential') && strcmp(domain.boundaries(4),'essential'))
%                 nc(4,i) = 0.5*(physical.Tb(3) + physical.Tb(4)); %We take the mean value of boundaries
%             elseif(strcmp(domain.boundaries(3),'essential'))
%                 nc(4,i) = physical.Tb(3);
%             elseif(strcmp(domain.boundaries(4),'essential'))
%                 nc(4,i) = physical.Tb(4);
%             else
%                 error('No essential boundary?');
%             end
%         elseif(strcmp(domain.boundaries(3),'natural') && strcmp(domain.boundaries(4), 'natural'))
%             %if both boundaries are natural, the node is natural
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = 0.5*(physical.qb(3) + physical.qb(4)); %we take the mean value
%         elseif(isnan(nc(3,i)))
%             error('Neither essential nor natural boundary');
%         end
%     end
    
%     if(i >= 1 && i < domain.nElX + 1) %lower boundary
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(1), 'essential'))
%             nc(3,i) = 0;
%             nc(4,i) = physical.Tb(1);
%         elseif(strcmp(domain.boundaries(1), 'natural'))
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = physical.qb(1);
%         else
%             error('Neither essential nor natural boundary');
%         end
%     end
%     
%     if(mod(i,domain.nElX + 1) == 0 && i ~= (domain.nElX + 1)*(domain.nElY + 1)) %right boundary
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(2), 'essential'))
%             nc(3,i) = 0;
%             nc(4,i) = physical.Tb(2);
%         elseif(strcmp(domain.boundaries(2), 'natural') && nc(3,i) ~= 0)
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = physical.qb(2);
%         else
%             error('Neither essential nor natural boundary');
%         end
%     end
%     
%     if(i > (domain.nElX + 1)*domain.nElY + 1 && i <= (domain.nElX + 1)*(domain.nElY + 1)) %upper boundary
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(3), 'essential'))
%             nc(3,i) = 0;
%             nc(4,i) = physical.Tb(3);
%         elseif(strcmp(domain.boundaries(3), 'natural') && nc(3,i) ~= 0)
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = physical.qb(3);
%         else
%             error('Neither essential nor natural boundary');
%         end
%     end
%     
%     if(mod(i, domain.nElX + 1) == 1 && i ~= 1) %left boundary
%         boundaryflag = true;
%         if(strcmp(domain.boundaries(4), 'essential'))
%             nc(3,i) = 0;
%             nc(4,i) = physical.Tb(4);
%         elseif(strcmp(domain.boundaries(4), 'natural') && nc(3,i) ~= 0)
%             %Assign equation number j
%             nc(3,i) = j;
%             j = j + 1;
%             nc(5,i) = physical.qb(4);
%         else
%             error('Neither essential nor natural boundary');
%         end
%     end
    
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
nc(4, physical.essentialNodes) = physical.Tb(physical.boundaryType);
nc(5, physical.naturalNodes) = physical.qb(~physical.boundaryType);

end

