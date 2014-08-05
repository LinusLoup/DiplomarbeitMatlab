function j_E = normal_flux(E_p,nodes,triangles,edges,edge_triangles,uS_values)
%NORMAL_FLUX berechnet den Normalenfluß über alle Kanten aus E_p mit
%Normalenvektor n. Hierbei ist E_p die Menge der Kanten, über die der
%Normalenfluss jeweils berechnet werden soll. Weiterhin sind nodes die
%Knoten der Triangulierung, triangles die Matrix mit den Dreiecksindizes,
%edges die Mittelpunkts-/Kantenindexmatrix, edge_triangles die
%Kanten-Dreicks-Zuordnungs-Matrix und uS_values die Funktionswerte der
%approximierten Lösung u_S an den Knoten.


%% Initialisierung:
j_E = zeros(size(E_p));
    
%% Berechnung der Nachbarn vom Kanten-Set E_p:
[neighbours,flag] = neighbourhood(E_p,edge_triangles,'edges');
        
%% Berechnung der einzelnen Normalenflüsse:
for k = 1:length(E_p)
               
    % Berechnung der Kantenpunkte:
    edge_poi_ind = edges(3:4,E_p(k));
    edge_poi = nodes(:,edge_poi_ind);
            
    % Berechnung der Gradienten von u_S auf T_1 und T_2:
    neigh_tri_ind = neighbours(:,k);
    neigh_tri = triangles(1:3,neigh_tri_ind);
    p_T = nodes(:,neigh_tri);
    uS_T = uS_values(neigh_tri);
            
    p_T1 = p_T(:,1:3);
    p_T2 = p_T(:,4:6);
    uS_T1 = uS_T(:,1);
            
    if flag(k) == 0
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
    
    normal = graduS_T*n;
    j_E(k) = normal(2)-normal(1);
end

end