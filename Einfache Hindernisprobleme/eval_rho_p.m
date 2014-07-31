function rho_p = eval_rho_p(nodes,triangles,midpoints,midtri,u_S,eps_V,fun)
%EVAL_RHO_P berechnet den lokalen Anteil rho_P von rho_S. Mitgegeben werden
%die Knoten (Punkte) nodes, die Dreiecke triangles, sowie die 
%Kantenmittelpunkte midpoints und die Kantenmittelpunkt-Dreiecks-Zuordnung
%midtri. Die Parameter eps_V und fun bezeichnen die Lösung des
%Defekt-Problems und die Funktion, die im Integral von rho_S enthalten ist.
%u_S sind die Funktionswerte der linearen Lösung des approximierten
%Problems.

% Initialisierung verwendeter Größen:
np = size(nodes,2);
ntri = size(triangles,2);
nmp = size(midpoints,2);
rho_p = zeros(np,1);
my_triangle = [triangles(1:3,:);zeros(1,ntri)];

% Hut- und Bubble-Funktionen:
phi_P = @(xi,eta) [1-xi-eta; xi; eta];
phi_E = @(xi,eta) [4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
% Berechnung der Gewichte und Funktionswerte für das Flächenintegral:
[wi,~,phi_E_values] = quad_tri([0,1,0;0,0,1],phi_E,7);
[~,~,phi_P_values] = quad_tri([0,1,0;0,0,1],phi_P,7);
        
% lokale Gewichte & Stützstellen für die Gaußquadratur für das Kurvenint.:
%wi_local = [1,1];
nodes_local = [-sqrt(1/3),sqrt(1/3)];
phiPE_local = @(xi) (xi+1)/2.*(-xi.^2+1); % (-xi+1)/2.*(-xi.^2+1)];
local_phiPE_int = sum(phiPE_local(nodes_local));

for k = 1:np
    tri_index = find(my_triangle==k);
    phi_p_local = mod(tri_index,4);
    w_p = ceil(tri_index/4);
    E_index = midtri(:,w_p);
    
    % Berechnung des Flächenintegrals von rho_p über w_p:
    for j = 1:length(w_p)
        % Punkte des jeweiligen Dreiecks, auf dem wir uns befinden, mit
        % Jacobi-Determinante:
        mypoi = nodes(:,triangles(1:3,w_p(j)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % Gaußpunkte und die lokalen Anteile von eps_V:
        eps_V_local = eps_V(E_index(:,j));
        [~,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % die Funktionswerte der jeweiligen lokalen Hutfunktion phi_P:
        phi_Pl_values = phi_P_values(phi_p_local(j),:);
        % die Berechnung des ersten Integrals von rho_p:
        rho_p(k) = rho_p(k) + J * sum(wi.*fun(gauss(1,:),gauss(2,:)).*...
            (eps_V_local'*phi_E_values).*(phi_Pl_values));
    end
    
    % Berechnung des Kurvenintegrals über die Kanten E\in E_p:
        % Berechnung der Menge von Kanten E_p auf denen der Punkt P liegt:
        my_midpoints = [midpoints(3:4,:);zeros(1,nmp)];
        help_Ep = find(my_midpoints==k);
        E_p = ceil(help_Ep/3);
    
        % Berechnung der Nachbarn vom Kanten-Set E_p:
        [neighbours,flag] = neighbourhood(E_p,midtri,'edges');
        %no_bound_edge = find(flag);
        %E_p = E_p(no_bound_edge);
        %neighbours = neighbours(:,no_bound_edge);
        
        % Berechnung des Integrals mittels Gaußquadratur:
        for i = 1:length(E_p)
            % Berechnung der Kantenpunkte:
            edge_poi_ind = midpoints(3:4,E_p(i));
            edge_poi = nodes(:,edge_poi_ind);
            
            % Berechnung der Gradienten von u_S auf T_1 und T_2:
            neigh_tri_ind = neighbours(:,i);
            neigh_tri = triangles(1:3,neigh_tri_ind);
            p_T = nodes(:,neigh_tri);
            uS_T = u_S(neigh_tri);
            
            p_T1 = p_T(:,1:3);
            p_T2 = p_T(:,4:6);
            uS_T1 = uS_T(:,1);
            
            if flag(i) == 0
                uS_T2 = zeros(3,1);
            else
                uS_T2 = uS_T(:,2);
            end
            
            graduS_T = [gradu(p_T1,uS_T1);gradu(p_T2,uS_T2)];
            
            % Berechnung des Normalenvektors zwischen T_1 und T_2:
            connect = edge_poi(:,2)-edge_poi(:,1);
            orth_connect = [-connect(2);connect(1)];
            n = 1/norm(orth_connect)*orth_connect;
            
                % Test, ob n von T_1 nach T_2 zeigt:
            p3_T1 = setdiff(p_T1',edge_poi','rows')';
            edge_test = (p3_T1-edge_poi(:,1))';
            
            if edge_test*n > 0
                n = -n;
            end
            
            % Der Normalenfluß j_E:
            j_E = normal_flux(graduS_T,n);
            
            % Transformation des lokalen Integrals auf die globale Kante
            % und Multiplikation mit der Funktionaldeterminante bzgl. des
            % Kurvenintegrals:
            laenge = norm(edge_poi(:,1)-edge_poi(:,2));
            global_phiPE_int = local_phiPE_int*1/2*laenge;
            
            % Bestimmung des Funktionswertes von eps_V(x_E):
            epsV_loc = eps_V(E_p(i));
            
            % Hinzufügen der Kurvenintegrale zur rho_p:
            rho_p(k) = rho_p(k)+j_E*epsV_loc*global_phiPE_int;
        end
end


end