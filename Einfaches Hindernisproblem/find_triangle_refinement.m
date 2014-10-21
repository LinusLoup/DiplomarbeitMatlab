function triangle_index = find_triangle_refinement(rho_p,rho_global,...
    osc_local,osc_global,triangles,theta)
%FIND_TRIANGLE_REFINEMENT Summary of this function goes here
%   Detailed explanation goes here

%% Initialisierung von den Indizes der berechneten Punkte bzw. Dreicksindizes:
np = length(rho_p);
point_index_rho = zeros(np,1);
point_index_osc = zeros(np,1);
rho_bound = zeros(np,1);
osc_bound = zeros(np,1);
upperbound = rho_global+osc_global;
triangle_index = [];
counter = 1;

%% Berechnung der Punkte, deren lokaler Fehleranteil zur Verfeinerung führt:
while 1
    new_rho = max(rho_p);
    new_points = find(abs(rho_p - new_rho)<0.001);
    rho_p(new_points) = 0;
    number_new_points = length(new_points);
    rho_bound(counter) = number_new_points * new_rho;%+...
        %sqrt(sum(osc_local(new_points)));
    point_index_rho(counter:counter+number_new_points-1) = new_points;
    counter = counter + number_new_points;
    
    % Abbruchkriterium für die Punktesuche:
    if sum(rho_bound) >= theta*rho_global
        break;
    end
end

counter = 1;

while 0
    new_osc = max(osc_local);
    new_point = find(abs(osc_local - new_osc)<0.001);
    osc_local(new_point) = 0;
    number_new_points = length(new_point);
    osc_bound(counter) = number_new_points * new_osc;
    point_index_osc(counter:counter+number_new_points-1) = new_point;
    counter = counter + number_new_points;
    
    % Abbruchkriterium für die Punktesuche:
    if sum(osc_bound) >= theta_osc*osc_global
        break;
    end
end

%% Berechnung der Dreiecksindizes für die Verfeinerung:
point_index_rho = setdiff(point_index_rho,0);
point_index_osc = setdiff(point_index_osc,0);
point_index = union(point_index_rho,point_index_osc);

for k = 1:length(point_index)
    help_index = neighbourhood(point_index(k),triangles,'point');
    triangle_index = union(triangle_index,help_index);
end

end