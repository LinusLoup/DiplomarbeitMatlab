function triangle_index = find_triangle_refinement(rho_p,rho_global,...
    triangles,theta)
%FIND_TRIANGLE_REFINEMENT Summary of this function goes here
%   Detailed explanation goes here

%% Initialisierung von den Indizes der berechneten Punkte bzw. Dreicksindizes:
np = length(rho_p);
point_index = zeros(np,1);
rho_local = zeros(np,1);
triangle_index = [];
%counter = 1;

%% Berechnung der Punkte, deren lokaler Fehleranteil zur Verfeinerung führt:
for i = 1:np
    [new_rho,new_point] = max(rho_p);
    %new_rho = max(rho_p);
    %new_point = find(rho_p == new_rho);
    rho_p(new_point) = 0;
    %number_new_points = length(new_point);
    %rho_local(counter) = number_new_points * new_rho;
    %point_index(counter:counter+number_new_points-1) = new_point;
    %counter = counter + number_new_points;
    rho_local(i) = new_rho;
    point_index(i) = new_point;
    
    % Abbruchkriterium für die Punktesuche:
    if sum(rho_local) >= theta*rho_global
        break;
    end
end

%% Berechnung der Dreiecksindizes für die Verfeinerung:
point_index = setdiff(point_index,0);

for k = 1:length(point_index)
    help_index = neighbourhood(point_index(k),triangles,'point');
    triangle_index = union(triangle_index,help_index);
end

end