function triangle_index = find_triangle_refinement(rho_p,rho_global,...
    triangles,theta)
%FIND_TRIANGLE_REFINEMENT Summary of this function goes here
%   Detailed explanation goes here

%% Initialisierung von den Indizes der berechneten Punkte bzw. Dreicksindizes:
np = length(rho_p);
point_index = zeros(np,1);
rho_local = zeros(np,1);
triangle_index = [];

%% Berechnung der Punkte, deren lokaler Fehleranteil zur Verfeinerung f�hrt:
for i = 1:np
    [new_rho,new_point] = max(rho_p);
    rho_p(new_point) = 0;
    rho_local(i) = new_rho;
    point_index(i) = new_point;
    
    % Abbruchkriterium f�r die Punktesuche:
    if sum(rho_local) >= theta*rho_global
        break;
    end
end

%% Berechnung der Dreiecksindizes f�r die Verfeinerung:
point_index = setdiff(point_index,0);

for k = 1:length(point_index)
    help_index = neighbourhood(point_index(k),triangles,'point');
    triangle_index = union(triangle_index,help_index);
end

end