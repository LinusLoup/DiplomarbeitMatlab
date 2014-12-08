function [J_error,rhoS_plot,IQ_plot,degree_of_freedom,time] = start_example4

clear
clear all 
clc

% loading of the geometry data (non-constant boundary):
data = load('mycylinder.mat');
% loadfunction f (here contant):
    function z = load1(x,y)
        z = zeros(size(x));
    end
vol_load_x = @(x,y) load1(x,y);

    function z = load2(x,y)
        z = zeros(size(x));
        index = find(y > 0);
        z(index) = -200*1/2*(1-y(index));
    end
vol_load_y = @(x,y) load2(x,y);

    function z = load3(x,y)
        z = zeros(size(y));
        
        index = find(x>-0.26 & x<0.26 & y>1.99);
        z(index) = -100*ones(size(index));
    end

surf_load_x = @(x,y) load1(x,y);
surf_load_y = @(x,y) load3(x,y);
% loading the exakt data for the given problem:
J_u = 0;
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
h = 0.05;
[p,e,t] = initmesh(data.mygeomg,'Hmax',h,'MesherVersion','R2013a');

% plot of the mesh with the nodes and midpoints/edges:
figure(1)
%subplot(2,1,1);
pdemesh(p,e,t);

% initialization of the global values:
u_S = [];
itermax = 2;          % maximum iteration depth
nmax = 1000000;          % maximum number of nodes
eps = 0.00001;           % upper bound for the hierarchical error estimate
theta_rho = 0.4;      % contraction parameter for local contributions of the error estimate

tic
% adaptive algorithm:
[u_S,p,e,t,midtri,midpoints,rhoS_plot,IQ_plot,J_error,recursion_depth, degree_of_freedom,time] = adaptive_refinement_solution(p,e,t,lambda,mu,u_S,vol_load_x, vol_load_y,surf_load_x,surf_load_y,my_obstacle, data,J_u,eps,theta_rho,nmax,itermax);
toc

% title('numbering of the nodes and triangles','FontSize',12);
% 
% ntri = size(t,2);
% np = size(p,2);
% 
% for j = 1 : ntri
%     tri = t(1:3,j);
%     poi = p(:,tri);
%     adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
%     text(adv(1),adv(2), num2str(j), 'FontSize',9);
% end
% 
% for i = 1 : np
%     text(p(1,i),p(2,i), num2str(i), 'FontSize',9);
% end
% 
% subplot(2,1,2);pdemesh(p,e,t);
% title('numbering of the edges/midpoints','FontSize',12);
%
% for j = 1 : ntri
%     tri = midtri(1:3,j);
%     poi = midpoints(:,tri);
%     adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
%     text(adv(1),adv(2), num2str(j), 'FontSize',9);
% end
% 
% for i = 1 : length(midpoints)
%     text(midpoints(1,i),midpoints(2,i), num2str(i), 'FontSize',9);
% end



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
pdeplot(p,e,t,'xydata',u_S(1,:),'zdata',u_S,'mesh','on','colormap','jet', 'colorbar','on');
%title('solution of the obstacle problem','FontSize',15)

figure(recursion_depth+4);
pdeplot(p,e,t,'xydata',u_S(2,:),'zdata',u_S,'mesh','on','colormap','jet', 'colorbar','on');
%title('solution of the obstacle problem','FontSize',15)

end