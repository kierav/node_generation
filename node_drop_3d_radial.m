function xyz = node_drop_3d_radial (R, ninit, dotmax, radius, vargin)

% --- Input parameters ---
%   R           Radius of sphere be filled by nodes; 
%   ninit       Number of PDP entries: [th_init,phi_init]
%   dotmax      Upper bound for number of nodes to place
%   radius      The function radius(xyz) provides exclusion radius to be used
%               at location (x,y,z).
%   vargin      Additional parameters for exclusion radius function
%
% --- Output parameter ---
%   xyz          Array xyz(:,3) with the generated node locations

dotnr   = 0;                            % Counter for the placed nodes
rng(0);                                 % Initialize random number generator
pdp     = 1e-4*rand(ninit);             % Array to hold PDPs
placednodes = zeros(dotmax,3);          % Array to store produced dot locations in spherical coordinates
nodeindices = zeros(size(pdp));         % Array of pointers to the produced node locations 
excess_height = 0.3;                    % Percentage of the height to go over
dtheta = (2*pi)/(ninit(1));             % Grid size
dphi = pi/(ninit(2)-1);
[rm,idr]  = max(pdp(:));                % Locate PDP with minimum r-coordinate
[i1,i2] = ind2sub(size(pdp),idr);

while rm <= R*(1+excess_height) && dotnr < dotmax  
    % --- Add new node to generated nodes
    dotnr = dotnr + 1;                  
    placednodes(dotnr,:) = [pdp(i1,i2),dtheta*(i1-1),-pi/2 + dphi*(i2-1),];    
    nodeindices(i1,i2) = dotnr;
    
    [x,y,z] = sph2cart(placednodes(dotnr,2),placednodes(dotnr,3),placednodes(dotnr,1));
    r = radius([x,y,z],vargin);    
    
    % --- Correction for minimal spacing between nodes 
%     [x2,y2,z2] = sph2cart(placednodes(dotnr,2),placednodes(dotnr,3),placednodes(dotnr,1)+r);
%     r2 = radius([x2,y2,z2],vargin);
%     r = max([r,r2]);
    
    % --- Find PDPs inside the new circle
    deltaphi = acos(1-r^2/2/pdp(i1,i2)^2);
    deltatheta = acos((1-r^2/2/pdp(i1,i2)^2-(sin(placednodes(dotnr,3)))^2)/(cos(placednodes(dotnr,3)))^2);
    
    ileft  = i1 - ceil(deltatheta/dtheta);    
    if ileft < 1
        ileft = ninit(1) + ileft;
    end
    iright = i1 + ceil(deltatheta/dtheta);
    if iright > ninit(1)
        iright = iright - ninit(1);
    end
    ibottom = i2 - ceil(deltaphi/dphi);
    if ibottom < 1
        % Cover all theta if at base of sphere
        ileft = 1;
        iright = ninit(1);
        ibottom = 1;
    end
    itop = i2 + ceil(deltaphi/dphi);
    if itop > ninit(2)
        % Cover all theta if at top of sphere
        ileft = 1;
        iright = ninit(1);
        itop = ninit(2);
    end
    
    if ileft < iright
        xx = ileft:iright;
    else
        xx = [1:iright,ileft:ninit(1)];
    end
    yy = ibottom:itop;
    [X,Y] = ndgrid(xx,yy);
    

    % --- Update heights of PDPs within radius 
    angleFactor = cos(placednodes(dotnr,3))*cos(-pi/2+dphi*(Y-1)).*cos(placednodes(dotnr,2)-dtheta*(X-1)) ...
                    + sin(placednodes(dotnr,3))*sin(-pi/2+dphi*(Y-1));
    height =  sqrt(placednodes(dotnr,1)^2*(angleFactor.^2-1)+r^2);
    pdp(xx,yy) = max(pdp(xx,yy),pdp(i1,i2)*angleFactor+real(height));
    
    % --- Identify next node location as a local minimum of the PDPs
    [rm,ix] = min(pdp(xx,yy));
    [rm,iy] = min(rm);
    
    i1 = ileft+ix(iy)-1;
    if i1 > ninit(1)
        i1 = i1 - ninit(1);
    end
    i2 = ibottom+iy-1;
    
    while i2 > 1 && i2 < ninit(2) 
        switch i1
            case 1 % Wrap around in theta
                [rm,ix] = min(pdp([ninit(1),1,2],i2-1:i2+1));
                [rm,iy] = min(rm);
                ix = ix(iy);
                if ix == 2 && iy == 2
                    break
                elseif ix == 1
                    i1 = ninit(1);
                    i2 = i2+iy-2;
                    rm = pdp(i1,i2);
                else
                    i1 = i1+ix-2;
                    i2 = i2+iy-2;
                    rm = pdp(i1,i2);
                end
            case ninit(1) % Wrap around in theta
                [rm,ix] = min(pdp([i1-1,i1,1],i2-1:i2+1));
                [rm,iy] = min(rm);
                ix = ix(iy);
                if ix == 2 && iy == 2
                    break
                elseif ix == 3
                    i1 = 1;
                    i2 = i2+iy-2;
                    rm = pdp(i1,i2);
                else
                    i1 = i1+ix-2;
                    i2 = i2+iy-2;
                    rm = pdp(i1,i2);
                end  
            otherwise
                [rm,ix] = min(pdp(i1-1:i1+1,i2-1:i2+1));
                [rm,iy] = min(rm);
                ix = ix(iy);
                if ix == 2 && iy == 2
                    break
                else
                    i1 = i1+ix-2;
                    i2 = i2+iy-2;
                    rm = pdp(i1,i2);
                end
        end
    end 
end                                  

placednodes = placednodes(1:dotnr,:);
placednodes = placednodes(placednodes(:,1)<=R,:);
[x,y,z] = sph2cart(placednodes(:,2),placednodes(:,3),placednodes(:,1));
xyz = [x y z];                   
