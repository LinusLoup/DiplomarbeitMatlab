function N0plus_set = N0plus(N0_set,obstacle,nodes,triangles,uS_values)
%N0PLUS berechnet die Menge N^{0+} zur Berechnung des Oszillationsterms
%osc_1.
%
%Hierbei sind N0_set ein Vektor der Punkteindizes für die Kontaktknoten,
%obstacle die Hindernisfunktion \psi, nodes die Knotenpunkte des Gitters,
%triangles die Dreiecke und uS_values ein Vektor aus den Funktionswerten
%der approximierten Lösung u_S.


%% Initialisierung von N0plus_set:
N0plus_set = zeros(size(N0_set));

%% Berechnung der Punkte auf Referenzdreieck zum Vergleich der Funktionen:
[xi,eta] = meshgrid(0:0.1:1,0:0.1:1);
index_out_of_bounds = find(xi+eta>1);
xi(index_out_of_bounds) = 0;
eta(index_out_of_bounds) = 0;

[m_mesh,n_mesh] = size(xi);
np = m_mesh * n_mesh;
index_within = setdiff(1:np,index_out_of_bounds);

%% Hut-Funktionen mit den Funktionswerten:
phi_P = @(xi,eta) [1-xi-eta; xi; eta];
phi_P_values = phi_P(xi,eta);

%% Berechnung der Knotenindizes von N0+:
for i = 1:length(N0_set)
    % Träger der Ansatzfunktion & Index der lokalen Ansatzfunktion:
    [phi_p_local,w_p] = find(triangles(1:3,:)==N0_set(i));
    flag = 0;                   % Zur Überprüfung der Bedingung der Menge
    
    for j = 1:length(w_p)
        % Eckpunkte des Dreiecks und Funktionswerte u_S:
        p_index = triangles(1:3,w_p(j));
        uS_local = uS_values(p_index);
        mypoi = nodes(:,p_index);
        x = mypoi(1,:);
        y = mypoi(2,:);
        
        % Trafo von lokalem auf globales Dreick: 
        xval = x(1)+(x(2)-x(1)).*xi+(x(3)-x(1)).*eta;
        yval = y(1)+(y(2)-y(1)).*xi+(y(3)-y(1)).*eta;
        
        % Berechnung der Funktionswerte:
        z = uS_local(1)*phi_P_values(1:m_mesh,:)+uS_local(2)*...
            phi_P_values(1+m_mesh:2*m_mesh,:)+uS_local(3)*...
            phi_P_values(1+2*m_mesh:3*m_mesh,:)-obstacle(xval,yval);
        
        % Überprüfung, ob u_S - phi > 0:
        switch phi_p_local(j)
            case 1
                new_within = setdiff(index_within,1);
            case 2
                new_within = setdiff(index_within,sub2ind([m_mesh,n_mesh],1,n_mesh));
            case 3
                new_within = setdiff(index_within,m_mesh);
        end
        
        %true_val = any(z(new_within) <= 0)
        if any(z(new_within) <= 0)
            flag = 1;
            break;
        end
    end
    
    if flag == 0
        N0plus_set(i) = N0_set(i);
    end
end

N0plus_set = setdiff(N0plus_set,0);

end