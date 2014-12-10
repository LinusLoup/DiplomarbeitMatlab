function [J_error,rhoS_plot,IQ_plot,degree_of_freedom,time] = start_example4

clear
clear all 
clc

% loading of the geometry data (non-constant boundary):
data = load('mycylinder.mat');
% loadfunction f (here contant):
    function z = zero_load(x,y)
        z = zeros(size(x));
    end

    % possible function for a volume load:
    function z = y_vol_load(x,y)
        z = zeros(size(x));
        index = find(y > 0);
        z(index) = -100*1/2*(1-y(index));
    end
    
    % function for the surface load:
    function z = y_surf_load(x,y)
        z = zeros(size(y));
        
        index = find(x>-0.26 & x<0.26 & abs(x.^2+(y-1).^2 - 1)<=0.001);
        z(index) = -100*ones(size(index));
    end

    % function for the pointwise node load:
    function z = node_load(x,y)
        z = zeros(size(y));
        
        index = find(x>-0.156434465040232 & x<0.156434465040231);
        z(index) = -100*ones(size(index))./length(index);
    end

vol_load_x = @(x,y) zero_load(x,y);
vol_load_y = @(x,y) zero_load(x,y);
surf_load_x = @(x,y) zero_load(x,y);
surf_load_y = @(x,y) zero_load(x,y);
point_load = @(x,y) node_load(x,y);
% loading the exakt data for the given problem:
J_u = -2.29;
% obstacle function:
    function z = my_obs(x,y)
        z = 1-sqrt(1-x.^2);
    end
my_obstacle = @(x,y) my_obs(x,y);

% mechanical values:
E = [7000,10^6];
nu = [0.3,0.45];
lambda = nu.*E./((1+nu).*(1-2*nu));
mu = E./(2.*(1+nu));

% initializing the mesh:
h = 0.17;
[p,e,t] = initmesh(data.mygeomg,'Hmax',h,'MesherVersion','R2013a');

% plot of the mesh with the nodes and midpoints/edges:
figure(1);
pdemesh(p,e,t);

% initialization of the global values:
u_S = [];
itermax = 7;           % maximum iteration depth
nmax = 30000;          % maximum number of nodes
eps = 0.0001;          % upper bound for the hierarchical error estimate
theta_rho = 0.35;      % contraction parameter for local contributions of the error estimate

tic
% adaptive algorithm:
[u_S,p,e,t,rhoS_plot,IQ_plot,J_error,recursion_depth, degree_of_freedom,time] = adaptive_refinement_solution(p,e,t,lambda,mu,u_S,vol_load_x, vol_load_y,surf_load_x,surf_load_y,point_load,my_obstacle, data,J_u,eps,theta_rho,nmax,itermax);
toc

% plot of the error and the error estimator:
figure(recursion_depth+2)
plot(1:recursion_depth,J_error,'--o',...
    1:recursion_depth,IQ_plot,'-.*',1:recursion_depth,rhoS_plot,'-.x');
ymin = min([min(J_error),min(IQ_plot),min(rhoS_plot)])-0.5;
ymax = max([max(J_error),max(IQ_plot),max(rhoS_plot)])+0.5;
axis([0.5,recursion_depth+0.5,ymin,ymax]);
legend('functional error','estimated error','error indicator','location',...
    'best');

% plot of the solution:
figure(recursion_depth+3);
pdeplot(p,e,t,'xydata',u_S(1,:),'zdata',u_S,'mesh','off','colormap', 'jet', 'colorbar','on');
title('solution (in x-direction) of the obstacle problem','FontSize',15)

figure(recursion_depth+4);
pdeplot(p,e,t,'xydata',u_S(2,:),'zdata',u_S,'mesh','off','colormap', 'jet', 'colorbar','on');
title('solution (in y-direction) of the obstacle problem','FontSize',15)

end