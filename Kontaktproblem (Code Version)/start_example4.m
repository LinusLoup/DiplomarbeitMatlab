function start_example4

clear
clear all 
clc

% loading of the geometry data (non-constant boundary):
data = load('mycylinder.mat');
% loadfunction f (here contant):
f_load_x = @(x,y) zeros(size(x));
f_load_y = @(x,y) -2*ones(size(x));
% loading the exakt data for the given problem:
% data_exact = load('fval_log_u.mat');
% J_u = data_exact.fval;
J_u = 0;
% obstacle function:
    function z = my_obs(x,y)
        z = zeros(size(y));
        
        index1 = find(y>=1-1/sqrt(2) & y<=1 & x<0);
        z(index1) = sqrt((x(index1)+1).^2+y(index1).^2);
        
        index2 = find(y>=1-1/sqrt(2) & y<=1 & x>=0);
        z(index2) = sqrt((x(index2)-1).^2+y(index2).^2);
        
        index3 = find(y<1-1/sqrt(2) & y>0);
        z(index3) = y(index3)./(1-y(index3));
        
        index4 = find(y==0);
        z(index4) = sqrt(x(index4).^2+1)-1;
    end
my_obstacle = @(x,y) my_obs(x,y);

% mechanical values:
E = 7000;
nu = 0.3;
lambda = nu*E/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

% initializing the mesh:
h = 0.5;
[p,e,t] = initmesh(data.mygeomg,'Hmax',h,'MesherVersion','R2013a');

% plot of the mesh with the nodes and midpoints/edges:
figure(1)
%subplot(2,1,1);
pdemesh(p,e,t);
title('numbering of the nodes and triangles','FontSize',12);

ntri = size(t,2);
np = size(p,2);

for j = 1 : ntri
    tri = t(1:3,j);
    poi = p(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : np
    text(p(1,i),p(2,i), num2str(i), 'FontSize',9);
end



% initialization of the global values:
u_S = [];
itermax = 6;          % maximum iteration depth
nmax = 1000000;          % maximum number of nodes
eps = 0.001;           % upper bound for the hierarchical error estimate
theta_rho = 0.4;      % contraction parameter for local contributions of the error estimate
theta_osc = 0.3;      % contraction parameter for local contributions of the oscillations

tic
% adaptive algorithm:
[u_S,p,e,t,midtri,midpoints,rhoS_plot,IQ_plot,J_error,osc_term, osc1_term,osc2_term,recursion_depth, degree_of_freedom] = adaptive_refinement_solution(p,e,t,lambda,nu,u_S,f_load_x,f_load_y,my_obstacle, data,J_u,eps,theta_rho,theta_osc,nmax,itermax);
toc

subplot(2,1,2);pdemesh(p,e,t);
title('numbering of the edges/midpoints','FontSize',12);

for j = 1 : ntri
    tri = midtri(1:3,j);
    poi = midpoints(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : length(midpoints)
    text(midpoints(1,i),midpoints(2,i), num2str(i), 'FontSize',9);
end

% plot of the error and the error estimator:
figure(2);
subplot(2,1,1);
plot(1:recursion_depth,osc1_term,':o',1:recursion_depth,osc2_term, '-.*',1:recursion_depth,osc_term,':x');
ymin = min([min(osc_term),min(osc1_term),min(osc2_term)])-5;
ymax = max([max(osc_term),max(osc1_term),max(osc2_term)])+5;
axis([0.5,recursion_depth+0.5,ymin,ymax]);
legend('osc1','osc2','oscillation','location','best');

subplot(2,1,2);
plot(1:recursion_depth,J_error,'--o',...
    1:recursion_depth,IQ_plot,'-.*',1:recursion_depth,rhoS_plot,'-.x');
ymin = min([min(J_error),min(IQ_plot),min(rhoS_plot)])-0.5;
ymax = max([max(J_error),max(IQ_plot),max(rhoS_plot)])+0.5;
axis([0.5,recursion_depth+0.5,ymin,ymax]);
legend('functional error','estimated error','error indicator','location',...
    'best');

% plot of the solution:
u_S = full(u_S);

figure(3);
pdeplot(p,e,t,'xydata',u_S,'zdata',u_S,'mesh','on','colormap','jet','colorbar','off');
%title('solution of the obstacle problem','FontSize',15)

figure(4);
pdeplot(p,e,t,'zdata',u_S);
%title('solution of the obstacle problem','FontSize',15)

end