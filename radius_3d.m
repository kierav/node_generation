function   r = radius_3d(xyz,rFactor) 
% Return exclusion radius at location (x,y,z)
% rFactor - scaling factor

r = (1.04 - exp(-((xyz(:,1)).^2+(xyz(:,2)).^2+(xyz(:,3)).^2)*3))/rFactor;
