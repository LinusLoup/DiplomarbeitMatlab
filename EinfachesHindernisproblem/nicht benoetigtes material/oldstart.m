clear 
clear all 
clc

%% Laden der Geometriedaten:
%data = load('mycircle.mat');   %Einheitskreis
data = load('square_with_unconst_dirichlet.mat');
%data = load('mysquare.mat');   %Quadrat: [-1,1]^2

%% Lastfunktion f:
%fun = @(x,y) zeros(1,length(x));

% konstante Lastfunktion:
% fun = @(x,y) -5*ones(size(x));
% data_exact = load('u_exact_const_f_new.mat');
% J_u = data_exact.fval;

% konstante Lastfunktion und nichtkonstanter Rand:
fun = @(x,y) -2*ones(size(x));
data_exact = load('fval_log_u.mat');
J_u = data_exact.fval;

% Parabelf�rmige Lastfunktion:
% fun = @(x,y) -18*x.^2-5*y.^2;
% data_exact = load('u_exact_parabel_f.mat');
% J_u = data_exact.fval;

%fun = @(x,y) -3*x.^2-5*y.^2;


%% Gitterinitialisierung:
h = 2;
%[p,e,t] = initmesh(data.mycircleg,'Hmax',h);   %Einheitskreis
[p,e,t] = initmesh(data.mysquareg,'Hmax',h);    %Quadrat: [-1,1]^2
[p,e,t] = refinemesh(data.mysquareg,p,e,t);

%% globale Initialisierungen:
refine_triangle = [];
u_S = [];
recursion_depth = 1;        % Rekursionstiefe
recmax = 20;                 % maximale Rekursionstiefe
nmax = 6000;                % maximale Anzahl der verwendeten Punkte
eps = 0.01;                 % obere Grenze f�r hierarchischen Fehlersch�tzer
theta_rho = 0.2;            % Schranke f�r lokalen und globalen Anteil vom FS
%theta_osc = 0.3;            % Schranke f�r lokalen zu globalem Anteil von Oszillation
rhoS_plot = zeros(recmax,1);% Vektor von rho_S in allen Rekursionsschritten
IQ_plot = zeros(recmax,1);  % Vektor mit dem hierarchischen Fehler -I_Q(eps_V)
J_error = zeros(recmax,1);  % Vektor mit Fehler zwischen den Funktionalen
osc_term = zeros(recmax,1); % Vektor mit den Oszillationstermen

nodes_vec = zeros(recmax,1);

tic

while 1
    %% Initialisierungen:
    ntri = size(t,2);
    np = size(p,2);
    
    nodes_vec(recursion_depth) = np;
    
    %% Berechnung der Mittelpunkte:
    [midpoints,midtri] = midpoints_of_triangle(t,p);
    nmp = size(midpoints,2);

    %% Berechnung der z-Werte vom Hindernis f�r Plot, sowie mit den Gitter- 
    %% und Mittelpunkten:
    %z_obs_prob = (-p(1,:).^2-p(2,:).^2+0.3)';
    %z_obs_midpoints = (-midpoints(1,:).^2-midpoints(2,:).^2+0.3)';
    %obstacle = @(x,y) -ones(size(x))';     % Hindernis z = -1
    obstacle = @(x,y) -zeros(size(x))';
    z_obs_prob = obstacle(p(1,:),p(2,:));
    z_obs_midpoints = obstacle(midpoints(1,:),midpoints(2,:));


    %% Assemblierung der Daten mit Berechnung der Dirichlet-Randdaten: 
    [A,f] = assemble(p,t,fun,7,'linear');
    %[A,f]=assempde(data.mycircleb,p,e,t,1,0,0);
    %[~,~,H,R]=assemb(data.mycircleb,p,e);
    [~,~,H,R]=assemb(data.mysquareb,p,e);


    %% L�sung der Variationsungleichung mit Active-Set- und Jacobi-Verfahren:
    % Active-Set-Methode:
    u_S = sparse(u_S);
    z_obs_prob = sparse(z_obs_prob);
    opts = optimset('Algorithm','interior-point-convex','LargeScale','on',...
       'Display','off');
    [u_S,J_uS] = quadprog(A,-f,[],[],H,R,z_obs_prob,[],u_S,opts);


% % Eliminieren von Dirichletpunkten aus der Matrix f�r das Jacobi-Verfahren:
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
    
    for j = 1 : nmp
        index = midpoints([3,4],j);
        u_S_mid(j) = (u_S(index(1))+u_S(index(2)))/2; 
    end
    
    %% L�sen des Defektproblems:
% Aufstellen der Dirichletrandbedingungs-Matrix H_Q: wird hier nicht
% ben�tigt, da die L�sung des Defektproblems am Rand ungleich 0 sein wird.
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
    %eps_V = quadprog(A_Q,-rhoS_phiE,[],[],[],[],z_obs_midpoints-u_S_mid,[]);
    
    % Probe durch exakte L�sung laut (2.10):
    [eps_V,rho_E,d_E,~] = exact_defect(p,t,midtri,rhoS_phiE,u_S_mid,...
        z_obs_midpoints);


% % Vergleich der ASM mit dem projektiven Jacobi-Verfahren:
% eps_V_jac = projected_jacobi(A_Q,rhoS_phiE,z_obs_midpoints-u_S_mid,...
%     rand(nmp,1),1e-15);
% 
% err_eps_V_quad = norm(eps_V-eps_V_exact);
% err_eps_V_jac = norm(eps_V_jac-eps_V_exact);
% fprintf('%s %s: %.20f.\n','Fehler zwischen exakter L�sung des',...
%     'Defektproblems und ASM-Lsg. ist',err_eps_V_quad);
% fprintf('%s %s: %.20f.\n','Fehler zwischen exakter L�sung des',...
%     'Defektproblems und Jacobi-Lsg. ist',err_eps_V_jac);


    %% Der hierarchische Fehlersch�tzer:
    % Berechnung von rho_s(eps_V) nach (3.5):
    rhoS_glob = eps_V'*rhoS_phiE;
    rhoS_plot(recursion_depth) = rhoS_glob;
    
    % Berechnung von -I_Q(eps_V) nach :
    IQ_plot(recursion_depth) = -1/2*(eps_V.^2)'*diag(A_Q) + rhoS_glob;
    
    % Berechnung des Fehlers J(u_S)-J(u) zwischen den Funktionalen:
    J_error(recursion_depth) = J_uS-J_u;

    % Berechnung des lokalen rho_p(eps_V) nach S.661 zwischen (3.9) und (3.10):
    rho_p = eval_rho_p(p,t,midpoints,midtri,u_S,eps_V,fun);
%    rho_E = eval_rho_E(p,t,midpoints,midtri,u_S,eps_V,fun);

    % Berechnen der Mengen N0, N0+, N+, N++, N0- f�r die Oszillationsterme:
    N0_set = N0(u_S,z_obs_prob);
    Nplus_set = Nplus(N0_set,p);
    [N0plus_set,N0minus_set] = N0plusminus(N0_set,p,t,midpoints,midtri,...
        fun,obstacle,u_S);
    N0minus_set;
    Nplusplus_set = Nplusplus(Nplus_set,midpoints,rho_E,d_E);
    
    % Berechnnung der Oszillationsterme:
    [osc1_term,osc1_local] = osc1(N0plus_set,z_obs_prob,p,t,u_S);
    [osc2_term,osc2_local] = osc2(Nplusplus_set,N0minus_set,p,t,midpoints,fun);
    osc_local = osc1_local + osc2_local;
    osc_term(recursion_depth) = osc1_term + osc2_term;

    % Bestimmung der zu verfeinernden Dreiecke:
    refine_triangle = find_triangle_refinement(rho_p,rhoS_glob,osc_local,...
        osc_term(recursion_depth),t,theta_rho);
   
    %% Verfeinerung des Gitters:
    [p_h,e_h,t_h,uS_h] = refinemesh(data.mysquareg,p,e,t,u_S,refine_triangle);
    
    %% Abbruchkriterium:
    % falls die Anzahl der Ecken �ber nmax liegt oder kein Dreieck mehr
    % verfeinert wird oder der hierarchische Fehlersch�tzer gen�gend klein 
    % ist.
    if (isempty(refine_triangle) || rhoS_glob < eps...
            || recursion_depth == recmax || length(p_h) > nmax)
        length(p_h)
        fprintf('%s %f.\n','Die Rekursionstiefe ist ',recursion_depth);
        break;
    else
        p = p_h;
        e = e_h;
        t = t_h;
        u_S = uS_h;
        recursion_depth = recursion_depth + 1;
    end
end

toc

%% Eliminieren der Nullen aus dem Vektoren der Fehler/-sch�tzer:
rhoS_plot = rhoS_plot(1:recursion_depth);
IQ_plot = IQ_plot(1:recursion_depth);
J_error = J_error(1:recursion_depth);
osc_term = osc_term(1:recursion_depth);


%% Plot vom Gitter, Eck- sowie Mittelpunkten:
figure(1)
subplot(2,2,1);pdemesh(p,e,t);
title('Nummerierung der Ecken und Dreiecke','FontSize',12);

for j = 1 : ntri
    tri = t(1:3,j);
    poi = p(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : np
    text(p(1,i),p(2,i), num2str(i), 'FontSize',9);
end


subplot(2,2,2);pdemesh(p,e,t);
title('Nummerierung der Kanten/Mittelpunkte','FontSize',12);

for j = 1 : ntri
    tri = midtri(1:3,j);
    poi = midpoints(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : length(midpoints)
    text(midpoints(1,i),midpoints(2,i), num2str(i), 'FontSize',9);
end

%% Plot der Fehler:
subplot(2,2,3:4);

plot(1:recursion_depth,J_error,'--o',1:recursion_depth,osc_term,':x',...
    1:recursion_depth,IQ_plot,'-.*');
ymin = min([min(J_error),min(IQ_plot),min(osc_term)])-10;
ymax = max([max(J_error),max(IQ_plot),max(osc_term)])+10;
axis([0.5,recursion_depth+0.5,ymin,ymax]);
legend('functional error','oscillations','estimated error','location',...
    'best');

%% Plot des Hindernisses:
figure(2);
[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
%z_obs = -x.^2-y.^2+0.3;
z_obs = obstacle(x,y);

subplot(2,1,1);surf(x,y,z_obs);
title('Hindernis','FontSize',15);


%% Plot der L�sungen:
u_S = full(u_S);
subplot(2,1,2); pdeplot(p,e,t,'zdata',u_S);
title('L�sung des Hindernisproblems','FontSize',15)


%% Ausgabe der Eigenwerte und Berechnung der Determinanten von A:
eig_val = eig(A);
%fprintf('%.25f\n',eig_val);
detA = det(A);

% Elimination von Rundungsfehlern in der Determinante von A:
if ~isempty(find(abs(eig_val) < 1e-10, 1))
    detA=0;
end

disp(['Determinante = ' num2str(detA)]);