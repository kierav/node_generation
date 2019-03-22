function  r = radius_trui(xy,~) 
% Return exclusion radius at location (x,y) based on trui image

persistent F

if isempty(F)                           % Read in the trui image once
    A  = double(imread('trui.png','PNG')); A = flipud(A(:,:,1)); 
    rf = @(s) 0.002+0.006*s+0.012*s.^8; % Conversion of brightness to exclusion radius 
    F  = rf(A/255);
end

ixy = round(255*xy);                    % Given a location (x,y), evaluate  
r = F(min(ixy(:,2),255)+256*min(ixy(:,1),255)+1); % the corresponding exclusion radius

