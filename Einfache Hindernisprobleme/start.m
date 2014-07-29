clear all

%% Laden der Geometriedaten:
%data = load('mycircle.mat');
data = load('mysquare.mat');

%% Lastfunktion f:
%fun = @(x,y) zeros(1,length(x));
fun = @(x,y) -5*ones(1,length(x));
%fun = @(x,y) 9*(-x.^2-y.^2);
%fun = @(x,y) -3*x-5*y.^2;


%% Gitterinitialisierung:
h = 1;
%[p,e,t] = initmesh(data.mycircleg,'Hmax',h);
%[p,e,t] = refinemesh(data.mycircleg,p,e,t);
%[p,e,t] = refinemesh(data.mycircleg,p,e,t);

[p,e,t] = initmesh(data.mysquareg,'Hmax',h);
% [p,e,t] = refinemesh(data.mysquareg,p,e,t);
% [p,e,t] = refinemesh(data.mysquareg,p,e,t);
% [p,e,t] = refinemesh(data.mysquareg,p,e,t);
triangle_flag = zeros(size(t,2),1);
u_S = [];
recursion_depth = 1;

tic

while 1
    %% Initialisierungen:
    [p,e,t,u_S] = refinemesh(data.mysquareg,p,e,t,u_S,find(triangle_flag));
    ntri = size(t,2);
    np = size(p,2);

    %% Berechnung und Plot der Mittelpunkte:
    [midpoints,midtri] = midpoints_of_triangle(t,p);
    nmp = size(midpoints,2);


    %% Berechnung der z-Werte vom Hindernis für Plot, sowie mit den Gitter- 
    %% und Mittelpunkten:
    %z_obs_prob = (-p(1,:).^2-p(2,:).^2+0.3)';
    %z_obs_midpoints = (-midpoints(1,:).^2-midpoints(2,:).^2+0.3)';
    obstacle = @(x,y) -ones(size(x))';
    z_obs_prob = obstacle(p(1,:),p(2,:));
    z_obs_midpoints = obstacle(midpoints(1,:),midpoints(2,:));


    %% Assemblierung der Daten mit Berechnung der Dirichlet-Randdaten: 
    [A,f] = assemble(p,t,fun,7,'linear');
    %[A,f]=assempde(data.mycircleb,p,e,t,1,0,0);
    %[~,~,H,R]=assemb(data.mycircleb,p,e);
    [~,~,H,R]=assemb(data.mysquareb,p,e);


    %% Lösung der Variationsungleichung mit Active-Set- und Jacobi-Verfahren:
    % Active-Set-Methode:
    [u_quadprog,fval] = quadprog(A,-f,[],[],H,R,z_obs_prob,[],u_S);


% % Eliminieren von Dirichletpunkten aus der Matrix für das Jacobi-Verfahren:
% non_bound_ind = find(prod(H == zeros(size(H))));
% ind_length = length(non_bound_ind);
% 
% A_dirichlet = A(non_bound_ind,non_bound_ind);
% f_dirichlet = f(non_bound_ind);
% 
% % projeziertes Jacobi-Verfahren:
% uu = projected_jacobi(A_dirichlet,f_dirichlet,z_obs_prob(non_bound_ind),...
%     rand(ind_length,1),1e-15);
% u_jacobi = dirichlet_boundary(uu,H,R);  %Einbauen der Dirichletrandbedingungen


    %% Berechnung der Funktionswerte von u_S auf den Mittelpunkten:
    u_S_mid = zeros(nmp,1);
    u_S = u_quadprog;
<<<<<<< HEAD

    
=======
  
>>>>>>> FETCH_HEAD
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
    [A_Q,rhoS_phiE] = assemble(p,t,fun,7,'bubble',u_S);
    eps_V = quadprog(A_Q,-rhoS_phiE,[],[],[],[],z_obs_midpoints-u_S_mid,[]);


% % Vergleich der ASM mit dem projektiven Jacobi-Verfahren:
% eps_V_jac = projected_jacobi(A_Q,rhoS_phiE,z_obs_midpoints-u_S_mid,...
%     rand(nmp,1),1e-15);
% 
% 
<<<<<<< HEAD
% % Probe durch exakte Lösung laut (2.10):
=======
 % Probe durch exakte Lösung laut (2.10):
>>>>>>> FETCH_HEAD
% [eps_V,a_phi] = exact_defect(p,t,midtri,rhoS_phiE,u_S_mid,z_obs_midpoints);
% 
% err_eps_V_quad = norm(eps_V-eps_V_exact);
% err_eps_V_jac = norm(eps_V_jac-eps_V_exact);
% fprintf('%s %s: %.20f.\n','Fehler zwischen exakter Lösung des',...
%     'Defektproblems und ASM-Lsg. ist',err_eps_V_quad);
% fprintf('%s %s: %.20f.\n','Fehler zwischen exakter Lösung des',...
%     'Defektproblems und Jacobi-Lsg. ist',err_eps_V_jac);


    %% Der hierarchische Fehlerschätzer:
    % Berechnung von rho_s(eps_V) nach (3.5):
    rhoS_glob = eps_V'*rhoS_phiE;

    % Berechnung des lokalen rho_p(eps_V) nach S.661 zwischen (3.9) und (3.10):
    rho_p = eval_rho_p(p,t,midpoints,midtri,u_S,eps_V,fun);

    break;
    
    % Bestimmung der zu verfeinernden Dreiecke:
<<<<<<< HEAD
    eps = 0.001;
=======
    eps = 0.0001;
>>>>>>> FETCH_HEAD
    theta = 0.005;
    triangle_flag = zeros(ntri,1);

    for k = 1:length(rho_p)
        if rho_p(k) >= theta*rhoS_glob
            neighbours = neighbourhood(k,t,'point');
            triangle_flag(neighbours) = 1;
        end
    end

<<<<<<< HEAD
=======
    % Abbruchskriterium: Falls kein Dreieck mehr verfeinert wird oder der
    % hierarchische Fehlerschätzer genügend klein ist.
>>>>>>> FETCH_HEAD
    if (sum(triangle_flag == zeros(ntri,1)) == ntri || rhoS_glob < eps)
        fprintf('%s %f.\n','Die Rekursionstiefe ist ',recursion_depth);
        break;
    end
    recursion_depth = recursion_depth + 1;
end

toc


%% Plot vom Gitter, Eck- sowie Mittelpunkten:
subplot(3,3,1);pdemesh(p,e,t);

for j = 1 : ntri
    tri = t(1:3,j);
    poi = p(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : np
    text(p(1,i),p(2,i), num2str(i), 'FontSize',9);
end


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

%% Plot des Hindernisses:
[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
%z_obs = -x.^2-y.^2+0.3;
z_obs = obstacle(x,y);

subplot(3,3,3);surf(x,y,z_obs);


%% Plot der Lösungen mit Active-Set-Methode & projiziertem Jacobi-Verfahren:
subplot(3,3,4:9); pdeplot(p,e,t,'zdata',u_S);
%subplot(3,3,[7,8,9]); pdeplot(p,e,t,'zdata',u_jacobi);
%subplot(3,3,[7,8,9]); pdesurf(p,t,u_quadprog);


% %% Berechnung des Fehlers zwischen Active-Set- und Jacobi-Verfahren:
% err_vec = (u_quadprog-u_jacobi);
% err = norm(err_vec);
% fprintf('%s %s: %.20f.\n','Fehler zwischen Active-Set-Methode und',...
%     'projiziertem Jacobi-Verfahren',err);
% %fprintf('Der Fehlervektor u_quadprog-u_jacobi ist:\n')
% %fprintf('%.10f\n.',err_vec);


%% Ausgabe der Eigenwerte und Berechnung der Determinanten von A:
eig_val = eig(A);
%fprintf('%.25f\n',eig_val);
detA = det(A);

% Elimination von Rundungsfehlern in der Determinante von A:
if ~isempty(find(abs(eig_val) < 1e-14, 1))
    detA=0;
end

disp(['Determinante = ' num2str(detA)]);