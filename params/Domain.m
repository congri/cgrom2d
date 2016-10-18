classdef Domain
    %class describing the finite element domain

    properties (SetAccess = private)
        
        nElX = 10;                  %number of finite elements in each direction
        nElY = 10;
        nEl                         %total number of elements
        nNodes                      %total number of nodes
        boundaryNodes               %Nodes on the domain boundary
        essentialNodes              %essential boundary nodes
        naturalNodes                %natural boundary nodes
        boundaryElements            %Elements on boundary, counterclockwise counted
        boundaryType                %true for essential, false for natural boundary node
        lx = 1;                     %domain size
        ly = 1;
        lElX                        %element size
        lElY
        AEl
        nEq                         %number of equations
        %lc gives node coordinates, taking in element number and local node
        %number
        lc
        %Nodal coordiante array
        %holds global nodal coordinates in the first two lines (x and y).
        %In the thrid line, the equation number is stored
        nodalCoordinates
        %globalNodeNumber holds the global node number, given the element number as row
        %and the local node number as column indices
        globalNodeNumber
        %total number of nodes
        totalNodeNumber
        %Shape function gradient array
        Bvec
        %essential boundary (yes or no) given local node and element number
        essentialBoundary
        %lm takes element number as row and local node number as column index
        %and gives equation number
        lm
        %Get mapping from equation number back to global node number
        id
        %equation number and local node number precomputation for sparse stiffness
        %assembly
        Equations
        LocalNode
        kIndex
  
    end
    
    methods
        
        function domainObj = Domain(nElX, nElY, lx, ly)
            %constructor
            if nargin > 0
                domainObj.nElX = nElX;
                if nargin > 1
                    domainObj.nElY = nElY;
                    if nargin > 2
                        domainObj.lx = lx;
                        if nargin == 4
                            domainObj.ly = ly;
                        end
                    end
                end
                if nargin > 4
                    error('Wrong number of input arguments')
                end
            end
            domainObj.lElX = domainObj.lx/domainObj.nElX;
            domainObj.lElY = domainObj.ly/domainObj.nElY;
            domainObj.nEl = domainObj.nElX*domainObj.nElY;
            domainObj.nNodes = (domainObj.nElX + 1)*(domainObj.nElY + 1);
            domainObj.AEl = domainObj.lElX*domainObj.lElY;
            domainObj.boundaryNodes = [1:(domainObj.nElX + 1),...
                2*(domainObj.nElX + 1):(domainObj.nElX + 1):(domainObj.nElX + 1)*(domainObj.nElY + 1),...
                ((domainObj.nElX + 1)*(domainObj.nElY + 1) - 1):(-1):((domainObj.nElX + 1)*domainObj.nElY + 1),...
                (domainObj.nElX + 1)*((domainObj.nElY - 1):(-1):1) + 1];
            domainObj.boundaryElements = [1:domainObj.nElX,...
                2*(domainObj.nElX):(domainObj.nElX):(domainObj.nElX*domainObj.nElY),...
                ((domainObj.nElX)*(domainObj.nElY) - 1):(-1):(domainObj.nElX*(domainObj.nElY - 1) + 1),...
                (domainObj.nElX)*((domainObj.nElY - 2):(-1):1) + 1];
            
            %local coordinate array. First index is element number, 2 is local node, 3 is x or y
            domainObj.lc = get_loc_coord(domainObj);
            domainObj.globalNodeNumber = get_glob(domainObj);
            domainObj.totalNodeNumber = domainObj.globalNodeNumber(end, end - 1);

        end
        
        function domainObj = setNodalCoordinates(domainObj, physical)
            domainObj.nodalCoordinates = get_coord(domainObj, physical);
            domainObj.lm = domainObj.globalNodeNumber;
            Sg = size(domainObj.globalNodeNumber);
            for i = 1:Sg(1)
                for j = 1:Sg(2)
                    domainObj.lm(i,j) = domainObj.nodalCoordinates(3,domainObj.globalNodeNumber(i,j));
                end
            end
            domainObj.id = get_id(domainObj.nodalCoordinates);
            [domainObj.Equations, domainObj.LocalNode] = get_equations(domainObj.nEl, domainObj.lm);
            domainObj.kIndex = sub2ind([4 4 domainObj.nEl], domainObj.LocalNode(:,1), domainObj.LocalNode(:,2), domainObj.LocalNode(:,3));
        end
        
        function domainObj = setBvec(domainObj)
            domainObj.nEq = max(domainObj.nodalCoordinates(3,:));
            %Gauss points
            xi1 = -1/sqrt(3);
            xi2 = 1/sqrt(3);
            
            domainObj.Bvec = zeros(8, 4, domainObj.nEl);
            for e = 1:domainObj.nEl
                for i = 1:4
                    domainObj.essentialBoundary(i,e) = ~isnan(domainObj.nodalCoordinates(4,domainObj.globalNodeNumber(e,i)));
                end
                %short hand notation
                x1 = domainObj.lc(e,1,1);
                x2 = domainObj.lc(e,2,1);
                y1 = domainObj.lc(e,1,2);
                y4 = domainObj.lc(e,4,2);
                
                %Coordinate transformation
                xI = 0.5*(x1 + x2) + 0.5*xi1*(x2 - x1);
                xII = 0.5*(x1 + x2) + 0.5*xi2*(x2 - x1);
                yI = 0.5*(y1 + y4) + 0.5*xi1*(y4 - y1);
                yII = 0.5*(y1 + y4) + 0.5*xi2*(y4 - y1);
                
                %Assuming bilinear shape functions here!!!
                B1 = [yI-y4 y4-yI yI-y1 y1-yI; xI-x2 x1-xI xI-x1 x2-xI];
                B2 = [yII-y4 y4-yII yII-y1 y1-yII; xII-x2 x1-xII xII-x1 x2-xII];
                %Do not forget cross terms
                B3 = [yI-y4 y4-yI yI-y1 y1-yI; xII-x2 x1-xII xII-x1 x2-xII];
                B4 = [yII-y4 y4-yII yII-y1 y1-yII; xI-x2 x1-xI xI-x1 x2-xI];
                
                %Note:in Gauss quadrature, the differential transforms as dx = (l_x/2) d xi. Hence
                %we take the additional factor of sqrt(A)/2 onto B
                domainObj.Bvec(:,:,e) = (1/(2*sqrt(domainObj.AEl)))*[B1; B2; B3; B4];
            end
        end
        
        function domainObj = setBoundaries(domainObj, natNodes)    
            %natNodes holds natural nodes counted counterclockwise around domain, starting in lower
            %left corner
            domainObj.boundaryType = true(1, 2*domainObj.nElX + 2*domainObj.nElY);
            domainObj.boundaryType(natNodes) = false;
            domainObj.essentialNodes = domainObj.boundaryNodes(domainObj.boundaryType);
            domainObj.naturalNodes = domainObj.boundaryNodes(~domainObj.boundaryType);
        end
    end
    
end
























