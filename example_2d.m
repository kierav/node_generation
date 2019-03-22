%% 2D node generation based on trui image
close all;

ninit = 5e4;    % recommended background grid is 0.1*minimum exclusion radius
dotmax = 5e5;

xy = node_drop_2d([0 1 0 1],ninit,dotmax,@radius_trui);

figure
plot(xy(:,1),xy(:,2),'.k')
axis square; axis equal