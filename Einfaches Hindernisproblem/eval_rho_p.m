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
rho_p = zeros(np,1);

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
    % Zuordnung Dreieck <-> Ref-Fkt & Träger von phi_P, d.h. Dreiecksindizes
    [phi_p_local,w_p] = find(triangles(1:3,:)==k);
    % Kantenindizes der jeweiligen Dreiecke
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
        [~,E_p] = find(midpoints(3:4,:)==k);
    
        % Berechnung der Normalenflüsse für alle Kanten aus E_p:
        j_E = normal_flux(E_p,nodes,triangles,midpoints,midtri,u_S);

        for i = 1:length(E_p)    
            % Berechnung der Kantenpunkte:
            edge_poi_ind = midpoints(3:4,E_p(i));
            edge_poi = nodes(:,edge_poi_ind);
            
            % Transformation des lokalen Integrals auf die globale Kante
            % und Multiplikation mit der Funktionaldeterminante bzgl. des
            % Kurvenintegrals:
            laenge = norm(edge_poi(:,1)-edge_poi(:,2));
            global_phiPE_int = local_phiPE_int*1/2*laenge;
            
            % Bestimmung des Funktionswertes von eps_V(x_E):
            epsV_loc = eps_V(E_p(i));
            
            % Hinzufügen der Kurvenintegrale zur rho_p:
            rho_p(k) = rho_p(k)+j_E(i)*epsV_loc*global_phiPE_int;
        end
end

end