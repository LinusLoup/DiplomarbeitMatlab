function [J_error,rhoS_plot,IQ_plot,degree_of_freedom,time] = start_example2

clear 
clear all 
clc

% loading of the geometry data (non-constant boundary):
data = load('mylshape.mat');
% loadfunction f (here contant):
    function z = my_fun(x,y)
        [m,n] = size(x);
        z = zeros(m,n);
        r = sqrt(x.^2+y.^2);
        r_new = 2*r-1/2;
        
        index_gamma1 = find(sqrt(x.^2+y.^2)>=1/4 & sqrt(x.^2+y.^2)<3/4 & y>=0);
        
        z(index_gamma1) = -(r(index_gamma1)).^(2/3)...
            .*sin(2/3*atan2(y(index_gamma1),x(index_gamma1)))...
            .*(1./r(index_gamma1).*(-60*r_new(index_gamma1).^4+120*...
            r_new(index_gamma1).^3-60*r_new(index_gamma1).^2)+(-480*...
            r_new(index_gamma1).^3+720*r_new(index_gamma1).^2-240*...
            r_new(index_gamma1)))-4/3*r(index_gamma1).^(1/3).*(-60*...
            r_new(index_gamma1).^4+120*r_new(index_gamma1).^3-60*...
            r_new(index_gamma1).^2).*sin(2/3*atan2(y(index_gamma1), x(index_gamma1)));
        
        index_gamma1 = find(sqrt(x.^2+y.^2)>=1/4 & sqrt(x.^2+y.^2)<3/4 & y<0);
        
        z(index_gamma1) = -(r(index_gamma1)).^(2/3)...
            .*sin(2/3*(atan2(y(index_gamma1),x(index_gamma1))+2*pi))...
            .*(1./r(index_gamma1).*(-60*r_new(index_gamma1).^4+120*...
            r_new(index_gamma1).^3-60*r_new(index_gamma1).^2)+(-480*...
            r_new(index_gamma1).^3+720*r_new(index_gamma1).^2-240*...
            r_new(index_gamma1)))-4/3*r(index_gamma1).^(1/3).*(-60*...
            r_new(index_gamma1).^4+120*r_new(index_gamma1).^3-60*...
            r_new(index_gamma1).^2).*sin(2/3*(atan2(y(index_gamma1), x(index_gamma1))+2*pi));
        
        index_gamma2 = find(sqrt(x.^2+y.^2)>5/4);
        z(index_gamma2) = -1;
    end
fun = @(x,y) my_fun(x,y);
% loading the exakt data for the given problem
J_u = -0.5;
% obstacle function:
my_obstacle = @(x,y) zeros(size(x));


% initializing the mesh:
h = 2;
[p,e,t]=initmesh(data.mygeomg,'Hmax',h,'jiggle','on','MesherVersion', 'R2013a');
[p,e,t] = refinemesh(data.mygeomg,p,e,t);
[p,e,t] = refinemesh(data.mygeomg,p,e,t);

% initialization of the global values:
u_S = [];
itermax = 3;          % maximum iteration depth
nmax = 25000;          % maximum number of nodes
eps = 0.001;           % upper bound for the hierarchical error estimate
theta_rho = 0.3;      % contraction parameter for local contributions of the error estimate
theta_osc = 0.3;      % contraction parameter for local contributions of the oscillations

tic
% adaptive algorithm:
[u_S,p,e,t,midtri,midpoints,rhoS_plot,IQ_plot,J_error,osc_term, osc1_term,osc2_term,recursion_depth,degree_of_freedom,time] = adaptive_refinement_solution(p,e,t,u_S,fun,my_obstacle, data,J_u,eps,theta_rho,theta_osc,nmax,itermax);
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
plot(1:recursion_depth,osc1_term,':o',1:recursion_depth,osc2_term, '-.*',1:recursion_depth,osc_term,':x');
ymin = min([min(osc_term),min(osc1_term),min(osc2_term)])-5;
ymax = max([max(osc_term),max(osc1_term),max(osc2_term)])+5;
axis([0.5,recursion_depth+0.5,ymin,ymax]);
legend('osc1','osc2','oscillation','location','best');

subplot(2,1,2);
plot(1:recursion_depth,J_error,'--o',...
    1:recursion_depth,IQ_plot,'-.*',1:recursion_depth,rhoS_plot,'-.x');
ymin = min([min(J_error),min(IQ_plot),min(rhoS_plot)])-0.02;
ymax = max([max(J_error),max(IQ_plot),max(rhoS_plot)])+0.02;
axis([0.5,recursion_depth+0.5,ymin,ymax]);
legend('functional error','estimated error','error indicator','location',...
    'best');

% plot of the solution:
u_S = full(u_S);

figure(3);
pdeplot(p,e,t,'xydata',u_S,'zdata',u_S,'mesh','on','colormap', 'jet','colorbar','off');
title('solution of the obstacle problem','FontSize',15)

figure(4);
pdeplot(p,e,t,'zdata',u_S);
title('solution of the obstacle problem','FontSize',15)

end
