function [A,f] = assemble(points,triangle,fun,num_of_nodes,option,u_S)
%ASSEMBLE evaluates the global matrix A and load vector f out of the given nodes, triangles, loadfunction, number of nodes and the Galerkin approximation u_S.

% ordering: triangles <-> midpoints:
[midpoints,midtriangle] = midpoints_of_triangle(triangle,points);

% initialising the dimension and the solution:
np = size(points,2);
nt = size(triangle,2);
nmp = size(midpoints,2);

if nargin == 4
    option = 'linear';
    u_S = zeros(np,1);
end

if nargin == 5
    u_S = zeros(np,1);
end

% beginning of the assembling:
switch lower(option)
    case {'linear'}
    % linear hatfunctionen on the reference-element:
    hat = @(xi,eta) [1-xi-eta; xi; eta];
    
    % initialising of the global values:
    A = sparse(np,np);
    f = sparse(np,1);
    my_tri = triangle(1:3,:);
    
    case {'bubble','quadratic'}
    % bubblefunctionen on the reference-element:
    hat = @(xi,eta) [4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
    
    % initialising of the global values:
    A = sparse(nmp,nmp);
    f = sparse(nmp,1);
    my_tri = midtriangle;
end

% loop over the triangles for the assembling:
for i = 1:nt
    poi = points(:,triangle(1:3,i));
    u_S_loc = u_S(triangle(1:3,i));
    tri = my_tri(:,i);
    
    % evaluation of the local linear stiffness matrix and the Jacobian:
    [S,fl,J] = local_mat(poi,u_S_loc,option);
        
    % Evaluating the local vector fl: (here is Gauss-quadrature over triangles used)
    [wi,gauss,ansatz_values] = quad_tri(poi,hat,num_of_nodes);
            
    % quadrature over a triangle:
    for k = 1:3
        fl(k) = fl(k)+J*sum(wi.*fun(gauss(1,:),gauss(2,:)).*...
            ansatz_values(k,:));
    end
    
    % assembling into the global stiffness matrix A and the load vector f:
    for k = 1:3
        for l = 1:3
            A(tri(k),tri(l)) = A(tri(k),tri(l))+S(k,l);
        end
    
        f(tri(k)) = f(tri(k))+fl(k);
    end
end

end

