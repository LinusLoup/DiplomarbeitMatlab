function [A,f] = assemble(points,edges,triangle,vol_load_x,vol_load_y,...
    surf_load_x,surf_load_y,node_load,neumann_points,lambda,mu,...
    num_of_nodes,option,u_S)
%ASSEMBLE evaluates the global matrix A and load vector f out of the given nodes, triangles, loadfunctions, Neumann-point indices, Lame-constants, number of nodes and the Galerkin approximation u_S.

% ordering: triangles <-> midpoints:
[midpoints,midtriangle] = midpoints_of_triangle(triangle,points);
% initialising the dimension and the solution:
np = size(points,2);
nt = size(triangle,2);
nmp = size(midpoints,2);

if nargin == 12
    option = 'linear';
    u_S = zeros(2*np,1);
end

if nargin == 13
    u_S = zeros(2*np,1);
end

hat_edge = @(xi) [1/2.*(1-xi); 1/2.*(1+xi)];

% beginning of the assembling:
switch lower(option)
    case {'linear'}
    % linear hatfunctionen on the reference-element:
    hat = @(xi,eta) [1-xi-eta; xi; eta];
    
    % initialising of the global values:
    A = sparse(2*np,2*np);
    f = sparse(2*np,1);
    my_tri = triangle(1:3,:);
    neumann_ind = neumann_points;
    
    case {'bubble','quadratic'}
    % bubblefunctionen on the reference-element:
    hat = @(xi,eta) [4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
    
    % initialising of the global values:
    A = sparse(2*nmp,2*nmp);
    f = sparse(2*nmp,1);
    my_tri = midtriangle;
    neumann_ind = find(prod(ismember(midpoints(3:4,:),neumann_points)));
end

% loop over the triangles for the assembling:
for i = 1:nt
    poi_ind = triangle(1:3,i);
    poi = points(:,poi_ind);
    u_S_loc_x = u_S(2*triangle(1:3,i)-1);
    u_S_loc_y = u_S(2*triangle(1:3,i));
    tri = my_tri(:,i);
    
    if any(poi(2,:)>0)
        my_mu = mu(1);
        my_lambda = lambda(1);
    else
        my_mu = mu(2);
        my_lambda = lambda(2);
    end
    
    % evaluation of the local linear stiffness matrix and the Jacobian:
    [S,fl,J] = local_mat(poi,my_lambda,my_mu,u_S_loc_x,u_S_loc_y,option);    
    
    % Evaluating the local vector fl: (here is Gauss-quadrature over triangles used)
    [wi,gauss,ansatz_values] = quad_tri(poi,hat,num_of_nodes);
            
    % quadrature over a triangle:
    for k = 1:3
        fl(2*k-1) = fl(2*k-1)+J*sum(wi.*vol_load_x(gauss(1,:),gauss(2,:)).*...
            ansatz_values(k,:));
        fl(2*k) = fl(2*k)+J*sum(wi.*vol_load_y(gauss(1,:),gauss(2,:)).*...
            ansatz_values(k,:));
    end
    
    local_neumann_ind = ismember(poi_ind,neumann_points);
    if sum(local_neumann_ind) == 2
        edge_poi_local_ind = find(local_neumann_ind);
        edge_poi_ind = poi_ind(edge_poi_local_ind);
        edge_poi = poi(:,edge_poi_local_ind);
        edge_ind = find(prod(ismember(edges(1:2,:), edge_poi_ind)));
        [wi_edge,gauss_edge,J_edge,ansatz_val_edge] = quad_edge(edge_poi,hat_edge,5);

        for k = 1:3
            if edges(1,edge_ind) == poi_ind(k)
                fl(2*k-1) = fl(2*k-1)+J_edge*sum(wi_edge.*...
                    surf_load_x(gauss_edge(1,:),gauss_edge (2,:))...
                    .*ansatz_val_edge(1,:));
                fl(2*k) = fl(2*k)+J_edge*sum(wi_edge.*...
                    surf_load_y(gauss_edge(1,:),gauss_edge (2,:))...
                    .*ansatz_val_edge(1,:));
            else
                fl(2*k-1) = fl(2*k-1)+J_edge*sum(wi_edge.*...
                    surf_load_x(gauss_edge(1,:),gauss_edge (2,:))...
                    .*ansatz_val_edge(2,:));
                fl(2*k) = fl(2*k)+J_edge*sum(wi_edge.*...
                    surf_load_y(gauss_edge(1,:),gauss_edge (2,:))...
                    .*ansatz_val_edge(2,:));
            end
        end
    end
    
    % assembling into the global stiffness matrix A and the load vector f:
    for k = 1:3
        for l = 1:3
            A(2*tri(k)-1:2*tri(k),2*tri(l)-1:2*tri(l)) = ...
                A(2*tri(k)-1:2*tri(k),2*tri(l)-1:2*tri(l)) + ...
                S(2*k-1:2*k,2*l-1:2*l);
        end
    
        f(2*tri(k)-1:2*tri(k)) = f(2*tri(k)-1:2*tri(k))+fl(2*k-1:2*k);
    end
    
end

switch lower(option)
    case {'bubble','quadratic'}
        points = midpoints;
end
f(2*neumann_ind) = f(2*neumann_ind) + node_load(points(1,neumann_ind),...
    points(2,neumann_ind))';
end