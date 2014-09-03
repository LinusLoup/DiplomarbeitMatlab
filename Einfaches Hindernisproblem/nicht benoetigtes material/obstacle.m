function [zgrid,z_pde_mesh,boundary] = obstacle(pde_mesh,geom)
    %% x-y-Gitter:
    [mx,my] = meshgrid(-1:0.1:1,-1:0.1:1);
    zgrid = zeros(size(mx));
     %% Distanz-Hindernis: (Pyramide)
     %% Randpunkte der Geometrie (Quadrat):
     points = boundary_points(geom);
    
%     %% Berechnung des Funktionswertes vom Hindernis: z = dist(x,dOmega)
%     for k = 1:length(mx)
%         for l = 1:length(my)
%             zgrid(k,l) = max(my_dist([mx(k,l);my(k,l)],points)-3/5,0);
%         end
%     end
     dist_of_mesh_to_bound = my_dist(pde_mesh,points);
    
    % Funktionswerte für das Netz der Triagulierung:
    z_pde_mesh = max(dist_of_mesh_to_bound-3/5,0);
    
    %% Indizes der Randpunkte:
    boundary = find(dist_of_mesh_to_bound < 10^(-3));

    %% Paraboloid:
    zgrid = max(-mx.^2-my.^2+0.5,0);
    z_pde_mesh = -pde_mesh(1,:).^2-pde_mesh(2,:).^2+0.3;

%     %% Kasten:
%     n = size(pde_mesh,2);
%     for i = 1:n
%         if dist_of_mesh_to_bound(i) >= 1/3
%             z_pde_mesh(i) = 2;
%         else
%             z_pde_mesh(i) = 0;
%         end
%     end
end