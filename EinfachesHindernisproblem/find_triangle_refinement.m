function triangle_index = find_triangle_refinement(rho_p,rho_global,osc_local,osc_global, triangles,theta_rho,theta_osc,option)
%FIND_TRIANGLE_REFINEMENT evaluates the possible indices of triangles, which will be refined, given the local contribution rho_p of the error indicator rho_global, the local contributions osc_local of the oscillation terms osc_global, the triangles and two parameter to provide a boundary.

% Initializing the indices of the evaluated points or triangleindices:
np = length(rho_p);
point_index = zeros(np,1);
bound = zeros(np,1);
triangle_index = [];
counter = 1;

if nargin == 7
    option = '';
end

switch lower(option)
    case 'symmetric'
        % evaluating the points, which have got a large contribution to the error indicator:
        while 1
            new_rho = max(rho_p);
            new_points = find(abs(rho_p - new_rho)<0.001);
            rho_p(new_points) = 0;
            number_new_points = length(new_points);
            bound(counter) = number_new_points * new_rho;
            point_index(counter:counter+number_new_points-1) = new_points;
            counter = counter + number_new_points;
    
            % termination criterion for the search of the points:
            if sum(bound) >= theta_rho*rho_global
                break;
            end
        end

        % evaluation of the points, which oscillation contribution is high:
        while 1
            % termination criterion for the search of the points:
            if sum(bound) >= theta_osc*osc_global
                break;
            end
            
            new_osc = max(osc_local);
            new_point = find(abs(osc_local - new_osc)<0.001);
            osc_local(new_point) = 0;
            number_new_points = length(new_point);
            bound(counter) = number_new_points * new_osc;
            point_index(counter:counter+number_new_points-1) = new_point;
            counter = counter + number_new_points;
        end
        
    otherwise
        % evaluating the points, which have got a large contribution to the error indicator:
        while 1
            [new_rho,index] = max(rho_p);
            rho_p(index) = 0;
            bound(counter) = new_rho;
            point_index(counter) = index;
            counter = counter + 1;
    
            % termination criterion for the search of the points:
            if sum(bound) >= theta_rho*rho_global
                break;
            end
        end

        % evaluation of the points, which oscillation contribution is high:
        while 1
            % termination criterion for the search of the points:
            if sum(bound) >= theta_osc*osc_global
                break;
            end
            
            [new_osc,index] = max(osc_local);
            osc_local(index) = 0;
            bound(counter) = new_osc;
            point_index(counter) = index;
            counter = counter + 1;
        end
end

% determination of the triangles indices for the refinement:
point_index = setdiff(point_index,0);

for k = 1:length(point_index)
    help_index = neighbourhood(point_index(k),triangles,'point');
    triangle_index = union(triangle_index,help_index);
end

end