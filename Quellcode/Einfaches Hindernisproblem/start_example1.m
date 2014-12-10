clear 
clear all 
clc

% loading of the geometry data (non-constant boundary):
data = load('square_with_unconst_dirichlet.mat');
% loadfunction f (here contant):
fun = @(x,y) -2*ones(size(x));
% loading the exakt data for the given problem:
data_exact = load('fval_log_u.mat');
J_u = data_exact.fval;
% obstacle function:
my_obstacle = @(x,y) zeros(size(x));

% initializing the mesh:
h = 2;
[p,e,t] = initmesh(data.mygeomg,'Hmax',h);    %square: [-1,1]^2

% initialization of the global values:
u_S = [];
itermax = 8;          % maximum iteration depth
nmax = 1000000;       % maximum number of nodes
eps = 0.001;          % upper bound for the hierarchical error estimate
theta_rho = 0.4;      % contraction parameter for local contributions of 
                      % the error estimate
theta_osc = 0.3;      % contraction parameter for local contributions of 
                      % the oscillations

tic
% adaptive algorithm:
[u_S,p,e,t,midtri,midpoints,rhoS_plot,IQ_plot,J_error,osc_term,...
    osc1_term,osc2_term,recursion_depth, degree_of_freedom] = ...
    adaptive_refinement_solution(p,e,t,u_S,fun,my_obstacle,data,J_u,eps,...
    theta_rho,theta_osc,nmax,itermax);
toc

% plot of the mesh with the nodes and midpoints/edges:
figure(1)
subplot(2,1,1);pdemesh(p,e,t);
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
plot(1:recursion_depth,osc1_term,':o',1:recursion_depth,osc2_term, '-.*',...
    1:recursion_depth,osc_term,':x');
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
pdeplot(p,e,t,'xydata',u_S,'zdata',u_S,'mesh','on','colormap',...
    'jet','colorbar','off');
title('solution of the obstacle problem','FontSize',15)

figure(4);
pdeplot(p,e,t,'zdata',u_S);
title('solution of the obstacle problem','FontSize',15)
