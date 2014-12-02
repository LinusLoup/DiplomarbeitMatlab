function [A,f] = assemble(points,triangle,load_x,load_y,lambda,mu,num_of_nodes,option,u_S)
%ASSEMBLE evaluates the global matrix A and load vector f out of the given nodes, triangles, loadfunction, number of nodes and the Galerkin approximation u_S.

% ordering: triangles <-> midpoints:
%[midpoints,midtriangle] = midpoints_of_triangle(triangle,points);

% initialising the dimension and the solution:
np = size(points,2);
nt = size(triangle,2);
%nmp = size(midpoints,2);

if nargin == 7
    option = 'linear';
    u_S = zeros(2*np,1);
end

if nargin == 8
    u_S = zeros(2*np,1);
end

% beginning of the assembling:
switch lower(option)
    case {'linear'}
    % linear hatfunctionen on the reference-element:
    hat = @(xi,eta) [1-xi-eta; xi; eta];
    
    % initialising of the global values:
    A = sparse(2*np,2*np);
    f = sparse(2*np,1);
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
    [S,fl,J] = local_mat(poi,lambda,mu,u_S_loc,option);
        
    % Evaluating the local vector fl: (here is Gauss-quadrature over triangles used)
    [wi,gauss,ansatz_values] = quad_tri(poi,hat,num_of_nodes);
            
    % quadrature over a triangle:
    for k = 1:3
        fl(2*k-1) = fl(2*k-1)+J*sum(wi.*load_x(gauss(1,:),gauss(2,:)).*...
            ansatz_values(k,:));
        fl(2*k) = fl(2*k)+J*sum(wi.*load_y(gauss(1,:),gauss(2,:)).*...
            ansatz_values(k,:));
    end
    
    % assembling into the global stiffness matrix A and the load vector f:
    for k = 1:3
        for l = 1:3
            A(2*tri(k)-1:2*tri(k),2*tri(l)-1:2*tri(l)) = A(2*tri(k)-1:2*tri(k),2*tri(l)-1:2*tri(l))+S(2*k-1:2*k,2*l-1:2*l);
        end
    
        f(2*tri(k)-1:2*tri(k)) = f(2*tri(k)-1:2*tri(k))+fl(2*k-1:2*k);
    end
end

end

