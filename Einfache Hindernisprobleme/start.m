%% Laden der Geometriedaten:
%data = load('mycircle.mat');
data = load('mysquare.mat');

%% Lastfunktion f:
%fun = @(x,y) zeros(1,length(x));
fun = @(x,y) -3*ones(1,length(x));
%fun = @(x,y) 5*(-x.^2-y.^2);
%fun = @(x,y) -3*x.^2-20*y.^2;


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
%[~,~,H,R]=assemb(data.mycircleb,p,e);
[~,~,H,R]=assemb(data.mysquareb,p,e);


%% Lösung der Variationsungleichung mit Active-Set- und Jacobi-Verfahren:
% Active-Set-Methode:
[u_quadprog,fval] = quadprog(A,-f,[],[],H,R,z_obs_prob,[]);


% Eliminieren von Dirichletpunkten aus der Matrix für das Jacobi-Verfahren:
non_bound_ind = find(prod(H == zeros(size(H))));
ind_length = length(non_bound_ind);

A_dirichlet = A(non_bound_ind,non_bound_ind);
f_dirichlet = f(non_bound_ind);

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
% Aufstellen der Dirichletrandbedingungs-Matrix H_Q: wird hier nicht
% benötigt, da die Lösung des Defektproblems am Rand ungleich 0 sein wird.
% % bound_ind = find(sum(H));   % --> Randindizes
% % l_bound_ind = length(bound_ind);
% % H_Q = zeros(l_bound_ind,nmp);
% % count_loop = 1;
% % 
% % for i = 1 : l_bound_ind
% %     for j = count_loop : nmp             % findet Randmittelpunkte
% %         if length(intersect(bound_ind,midpoints([3,4],j))) == 2
% %             H_Q(i,j) = 1;
% %             count_loop = j+1;
% %             break;
% %         end
% %     end
% % end

% Assemblierung der Matrix mit den Bubble-Fktn. und ASM:
[A_Q,rho_s] = assemble(p,t,fun,7,'bubble',u_S);
eps_V = quadprog(A_Q,-rho_s,[],[],[],[],z_obs_midpoints-u_S_mid,[]);

toc


% Vergleich der ASM mit dem projektiven Jacobi-Verfahren:
eps_V_jac = projected_jacobi(A_Q,rho_s,z_obs_midpoints-u_S_mid,...
    rand(nmp,1),1e-15);


% Probe durch exakte Lösung laut (2.10):
[eps_V_exact,a_phi] = exact_defect(p,t,midtri,rho_s,u_S_mid,z_obs_midpoints);

err_eps_V_quad = norm(eps_V-eps_V_exact);
err_eps_V_jac = norm(eps_V_jac-eps_V_exact);
fprintf('%s %s: %.20f\n.','Fehler zwischen exakter Lösung des',...
    'Defektproblems und ASM-Lsg. ist',err_eps_V_quad);
fprintf('%s %s: %.20f\n.','Fehler zwischen exakter Lösung des',...
    'Defektproblems und Jacobi-Lsg. ist',err_eps_V_jac);


%% Der hierarchische Fehlerschätzer:
% Berechnung von rho_s(eps_V) nach (3.5):
rhoS_glob = eps_V'*rho_s;

% Berechnung des lokalen rho_p(eps_V) nach S.661 zwischen (3.9) und (3.10):
rho_p = zeros(np,1);
my_triangle = [t(1:3,:);zeros(1,ntri)];
phi_P = @(xi,eta) [1-xi-eta; xi; eta];
phi_E = @(xi,eta) [4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
[wi,~,phi_E_values] = quad_tri([0,1,0;0,0,1],phi_E,7);
[~,~,phi_P_values] = quad_tri([0,1,0;0,0,1],phi_P,7);
        
for k = 1:np
    tri_index = find(my_triangle==k);
    phi_p_local = mod(tri_index,4);
    w_p = ceil(tri_index/4);
    E_index = midtri(:,w_p);
    
    % Berechnung des Flächenintegrals von rho_p über w_p:
    for l = 1:length(w_p)
        % Punkte des jeweiligen Dreiecks, auf dem wir uns befinden, mit
        % Jacobi-Determinante:
        mypoi = p(:,t(1:3,w_p(l)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % Gaußpunkte und die lokalen Anteile von eps_V:
        eps_V_local = eps_V(E_index(:,l));
        [~,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % die Funktionswerte der jeweiligen lokalen Hutfunktion phi_P:
        phi_Pl_values = phi_P_values(phi_p_local(l),:);
        % die Berechnung des ersten Integrals von rho_p:
        rho_p(k) = rho_p(k) + J * sum(fun(gauss(1,:),gauss(2,:)).*...
            (eps_V_local'*phi_E_values).*phi_Pl_values);
    end
    
    % Berechnung des Kurvenintegrals über die Kanten E\in E_p:
        % Berechnung der Menge von Kanten E_p auf denen der Punkt P liegt:
        my_midpoints = [midpoints(3:4,:);zeros(1,nmp)];
        help_Ep = find(my_midpoints==k);
        E_p = ceil(help_Ep/3);
    
        % Berechnung der Nachbarn vom Kanten-Set E_p:
        [old_neighbours,flag] = edge_neighbourhood(E_p,midtri);
        no_bound_edge = find(flag);
        E_p = E_p(no_bound_edge);
        neighbours = old_neighbours(:,no_bound_edge);
end


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