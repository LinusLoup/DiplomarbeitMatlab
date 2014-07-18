%% Laden der Geometriedaten:
%data = load('mycircle.mat');
data = load('mysquare.mat');

%% Lastfunktion f:
%fun = @(x,y) zeros(1,length(x));
fun = @(x,y) -3*ones(1,length(x));
%fun = @(x,y) 5*(-x.^2-y.^2);
%fun = @(x,y) -3*x.^2-10*y.^2;


%% Gitterinitialisierung mit Plot:
h = 2;
%[p,e,t] = initmesh(data.mycircleg,'Hmax',h);
%[p,e,t] = refinemesh(data.mycircleg,p,e,t);
%[p,e,t] = refinemesh(data.mycircleg,p,e,t);

[p,e,t] = initmesh(data.mysquareg,'Hmax',h);
%[p,e,t] = refinemesh(data.mysquareg,p,e,t);
%[p,e,t] = refinemesh(data.mysquareg,p,e,t);
subplot(3,3,1);pdemesh(p,e,t);
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

%% Berechnung und Plot der Mittelpunkte:
[midpoints,midtri] = midpoints_of_triangle(t,p);
nmp = size(midpoints,2);
subplot(3,3,2);pdemesh(p,e,t);

for j = 1 : ntri
    tri = midtri(1:3,j);
    poi = midpoints(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : length(midpoints)
    text(midpoints(1,i),midpoints(2,i), num2str(i), 'FontSize',9);
end


%% Berechnung der z-Werte vom Hindernis für Plot, sowie mit den Gitter- 
%% und Mittelpunkten:
[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
z_obs = -x.^2-y.^2+0.3;
subplot(3,3,3);surf(x,y,z_obs);

z_obs_prob = (-p(1,:).^2-p(2,:).^2+0.3)';
z_obs_midpoints = (-midpoints(1,:).^2-midpoints(2,:).^2+0.3)';


tic

%% Assemblierung der Daten mit Berechnung der Dirichlet-Randdaten: 
[A,f] = assemble(p,t,fun,7,'linear');
%[A,f]=assempde(data.mycircleb,p,e,t,1,0,0);
%[Q,G,H,R]=assemb(data.mycircleb,p,e);
[~,~,H,R]=assemb(data.mysquareb,p,e);


%% Eliminieren von Dirichletpunkten aus der Matrix für das Jacobi-Verfahren:
non_bound_ind = find(prod(H == zeros(size(H))));
ind_length = length(non_bound_ind);


A_dirichlet = A(non_bound_ind,non_bound_ind);
f_dirichlet = f(non_bound_ind);

%% Lösung der Variationsungleichung mit Active-Set- und Jacobi-Verfahren:
% Active-Set-Methode:
[u_quadprog,~] = quadprog(A,-f,[],[],H,R,z_obs_prob,[]);

% projeziertes Jacobi-Verfahren:
uu = projected_jacobi(A_dirichlet,f_dirichlet,z_obs_prob(non_bound_ind),...
    rand(ind_length,1),1e-15);
u_jacobi = dirichlet_boundary(uu,H,R);  %Einbauen der Dirichletrandbedingungen


%% Berechnung der Funktionswerte von u_S auf den Mittelpunkten:
u_S_mid = zeros(nmp,1);
u_S = u_quadprog;

for j = 1 : nmp
   index = midpoints([3,4],j);
   u_S_mid(j) = (u_S(index(1))+u_S(index(2)))/2; 
end

%% Lösen des Defektproblems:
% Aufstellen der Dirichletrandbedingungs-Matrix H_Q:
bound_ind = find(sum(H));   % --> Randindizes
l_bound_ind = length(bound_ind);
H_Q = zeros(l_bound_ind,nmp);
count_loop = 1;

for i = 1 : l_bound_ind
    for j = count_loop : nmp             % findet Randmittelpunkte
        if length(intersect(bound_ind,midpoints([3,4],j))) == 2
            H_Q(i,j) = 1;
            count_loop = j+1;
            break;
        end
    end
end

% Assemblierung der Matrix mit den Bubble-Fktn. und ASM:
[A_Q,rho_s] = assemble(p,t,fun,7,'bubble',u_S);
eps_V = quadprog(A_Q,-rho_s,[],[],H_Q,R,z_obs_midpoints-u_S_mid,[]);

toc


% Vergleich der ASM mit dem projektiven Jacobi-Verfahren:
non_bound_ind = find(prod(H_Q == zeros(size(H_Q))));
ind_length = length(non_bound_ind);

A_dirichlet_Q = A_Q(non_bound_ind,non_bound_ind);
f_dirichlet_Q = rho_s(non_bound_ind);

eps_V_jac = projected_jacobi(A_dirichlet_Q,f_dirichlet_Q,...
    z_obs_midpoints(non_bound_ind)-u_S_mid(non_bound_ind),...
    rand(ind_length,1),1e-15);
eps_V_jac = dirichlet_boundary(eps_V_jac,H_Q,R);


% Probe durch exakte Lösung:



%% Plot der Lösungen mit Active-Set-Methode & projiziertem Jacobi-Verfahren:
subplot(3,3,[4,5,6]); pdeplot(p,e,t,'zdata',u_quadprog);
%subplot(3,3,[7,8,9]); pdeplot(p,e,t,'zdata',u_jacobi);
subplot(3,3,[7,8,9]); pdesurf(p,t,u_quadprog);


%% Berechnung des Fehlers zwischen Active-Set- und Jacobi-Verfahren:
err_vec = (u_quadprog-u_jacobi);
err = norm(err_vec);
fprintf('%s %s: %.20f\n.','Fehler zwischen Active-Set-Methode und',...
    'projiziertem Jacobi-Verfahren',err);
%fprintf('Der Fehlervektor u_quadprog-u_jacobi ist:\n')
%fprintf('%.10f\n.',err_vec);


%% Ausgabe der Eigenwerte und Berechnung der Determinanten von A:
eig_val = eig(A);
%fprintf('%.25f\n',eig_val);
detA = det(A);

% Elimination von Rundungsfehlern in der Determinante von A:
if ~isempty(find(abs(eig_val) < 1e-14, 1))
    detA=0;
end

disp(['Determinante = ' num2str(detA)]);










%% Analytische Lösung (?) ohne Last (f=0) auf einem Kreis:
% 
% a = 0.3;
% b = 1;
% xhat = 1-sqrt(1-a/b);
% alpha = (a-b*xhat^2)/(1-xhat);
% u_ana = zeros(np,1);
% for i = 1:np
%     r = sqrt(p(1,i)^2+p(2,i)^2);
%    if r >= -xhat && r <= xhat
%        u_ana(i) = a-b*r^2;
%    elseif r < -xhat
%        u_ana(i) = alpha*(1+r);
%    else
%        u_ana(i) = alpha*(1-r);
%    end
% end
% 
% subplot(3,2,[5,6]); pdeplot(p,e,t,'zdata',u_ana);
% 
% fval
%1/2*u_quadprog'*A*u_quadprog - f'*u_quadprog
%1/2*u_jacobi'*A*u_jacobi - f'*u_jacobi
% 1/2*u_ana'*A*u_ana