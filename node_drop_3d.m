function xyz = node_drop_3d (box, ninit, dotmax, radius, vargin)
% NODE_DROP_3D generates quasi-uniform nodes in a 3-D bounded box with
% spacing between nodes specified by an exclusion radius function
%
% Copyright (C) 2019 Kiera van der Sande
%
% --- Input parameters ---
%   box         Size of box be filled by nodes: [xmin,xmax,ymin,ymax,zmin,zmax]
%   ninit       Number of PDP entries: [xinit,yinit]
%   dotmax      Upper bound for number of nodes to place
%   radius      The function radius(xyz) provides exclusion radius to be used
%               at location (x,y,z).
%   vargin      Additional parameters for exclusion radius function
%
% --- Output parameter ---
%   xyz         Array xyz(:,3) with the generated node locations

dotnr   = 0;                            % Counter for the placed nodes
rng(0);                                 % Initialize random number generator
xyz      = zeros(dotmax,3);             % Array to store produced node locations
nodeindices = zeros(size(pdp));         % Array of pointers to the produced node locations 
excess_height = 0.1;                    % Percentage of the height to go over
dx = (box(2)-box(1))/(ninit(1)-1);      % Grid size
dy = (box(4)-box(3))/(ninit(2)-1);
xx = box(1):10*dx:box(2); yy = box(3):10*dy:box(4);
[XX,YY] = meshgrid(xx,yy);
r = radius([XX(:),YY(:),box(1)*ones(size(XX(:)))],box(2),rFunction,rFactor);
pdp = box(5)+0.01*min(r)*rand(ninit);   % Array to hold PDPs
[zm,idx]  = min(pdp(:));                % Locate PDP with lowest z-coordinate
[i1,i2] = ind2sub(size(pdp),idx);

while zm <= (1+excess_height)*box(6) && dotnr < dotmax  
    % --- Add new node to generated nodes
    dotnr = dotnr + 1;                  
    xyz(dotnr,:) = [box(1)+dx*(i1-1),box(3)+dy*(i2-1),pdp(i1,i2)];    
    nodeindices(i1,i2) = dotnr;
    
    r = radius(xyz(dotnr,:),vargin);            
    
    % --- Find PDPs inside the new circle
    ileft  = max(1,i1 - floor(r/dx));
    iright = min(ninit(1),i1 + floor(r/dx));
    ibottom = max(1,i2 - floor(r/dy));
    itop = min(ninit(2),i2 + floor(r/dy));
    
    xx = ileft:iright;
    yy = ibottom:itop;
    [X,Y] = ndgrid(xx,yy);
    
    % --- Update heights of PDPs within radius 
    height =  sqrt(r^2 - (dx*(X-i1)).^2 - (dy*(Y-i2)).^2);
    pdp(ileft:iright,ibottom:itop) = (max(pdp(ileft:iright,ibottom:itop),pdp(i1,i2) + real(height)));
    
    % --- Identify next node location as a local minimum of the PDPs
    [zm,ix] = min(pdp(ileft:iright,ibottom:itop));
    [zm,iy] = min(zm);
    i1 = ileft+ix(iy)-1;
    i2 = ibottom+iy-1;
    
    searchr = min(2*ceil(r/dx),floor(ninit(1)/2)-1);
    
    while 1
        % Wrap around if a boundary is reached      
        if i1-searchr < 1
            xsearch = [ninit(1)+i1-searchr:ninit(1),1:i1+searchr];
        elseif i1+searchr > ninit(1)
            xsearch = [i1-searchr:ninit(1),1:i1+searchr-ninit(1)];
        else
            xsearch = i1-searchr:i1+searchr;
        end
        if i2-searchr < 1
            ysearch = [ninit(2)+i2-searchr:ninit(2),1:i2+searchr];
        elseif i2+searchr > ninit(2)
            ysearch = [i2-searchr:ninit(2),1:i2+searchr-ninit(2)];
        else
            ysearch = i2-searchr:i2+searchr;
        end
        [zm,ix] = min(pdp(xsearch,ysearch));
        [~,iy] = min(zm);
        ix = ix(iy);
        i1 = xsearch(ix);
        i2 = ysearch(iy);
        zm = pdp(i1,i2);
        
        % Stop once a local min has been found within the search radius
        if ix > searchr/2 && ix < length(xsearch)-searchr/2 && ...
                iy > searchr/2 && iy < length(ysearch)-searchr/2
            break
        end
    end  
end                                     

xyz = xyz(1:dotnr,:);                      % Remove unused entries in array xy
xyz = xyz(xyz(:,3)<=box(6),:);             % Remove any nodes placed above the box