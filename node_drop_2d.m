function xy = node_drop_2d (box, ninit, dotmax, radius, vargin)
% NODE_DROP_2D generates quasi-uniform nodes in a 2-D bounded box with 
% spacing between nodes specified by an exclusion radius function
%
% Copyright (C) 2019 Kiera van der Sande
%
% --- Input parameters ---
%   box         Size of box be filled by nodes; [xmin,xmax,ymin,ymax]
%   ninit       Upper limit on number of PDP entries
%   dotmax      Upper bound for number of nodes to place
%   radius      The function radius(xy) provides exclusion radius to be used
%               at location (x,y).
%   vargin      Additional parameters for exclusion radius function
%
% --- Output parameter ---
%   xy          Array xy(:,2) with the generated node locations

if nargin < 5
    vargin = [];
end

dotnr   = 0;                            % Counter for the placed nodes
rng(0);                                 % Initialize random number generator
pdp     = [linspace(box(1),box(2),ninit)',box(3)+ 1e-4*(box(4)-box(3))*rand(ninit,1)]; % Array to hold PDPs
xy      = zeros(dotmax,2);              % Array to store produced node locations
nodeindices = zeros(length(pdp),1);     % Array of pointers to the produced node locations 
excessheight = 0.1;                     % Percentage of the height to go over
dx      = pdp(2,1)-pdp(1,1);            % Grid size
[ym,i]  = min(pdp(:,2));                % Locate PDP with lowest y-coordinate

while ym <= (1+excessheight)*box(4) && dotnr < dotmax 
    % --- Add new node to generated nodes
    dotnr = dotnr + 1;                  
    xy(dotnr,:) = pdp(i,:);             
    nodeindices(i) = dotnr;             
    
    r = radius(xy(dotnr,:),vargin);     
    
    % --- Find PDPs inside the new circle
    ileft  = i - floor(r/dx);
    ileft = max(1,ileft);
    iright = i + floor(r/dx);
    iright = min(ninit,iright);

    % --- Update heights of PDPs within radius 
    pdp(ileft:iright,2) = max([pdp(ileft:iright,2), sqrt(r^2-(pdp(ileft:iright,1)-pdp(i,1)).^2)+pdp(i,2)],[],2);
    
    % --- Identify next node location as a local minimum of the PDPs
    if pdp(ileft,2)<pdp(iright,2)
        i = ileft;
    else
        i = iright;
    end
    
    searchr = min(2*ceil(r/dx),floor(ninit/2)-1);     
    
    while 1
        % Wrap around if a boundary is reached
        if i-searchr < 1 
            xsearch = [ninit+i-searchr:ninit,1:i+searchr];
        elseif i+searchr > ninit
            xsearch = [i-searchr:ninit,1:i+searchr-ninit];
        else
            xsearch = i-searchr:i+searchr;
        end
        [ym,ix] = min(pdp(xsearch,2));
        i = xsearch(ix);
        
        % Stop once a local min has been found within the search radius
        if ix > searchr/2 && ix < length(xsearch)-searchr/2
            break
        end
    end
end                                     

xy = xy(1:dotnr,:);                     % Remove unused entries in array xy
xy = xy(xy(:,2)<=box(4),:);             % Remove any nodes placed above the box